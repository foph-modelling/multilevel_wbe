library(patchwork)
library(tidyverse)
library(data.table)
library(dtwclust)
library(igraph)
model = readRDS(paste0("../",controls$savepoint,"ma5.3.2.rds"))

ww1 = readRDS(file=paste0("../",controls$savepoint,"ww_all.rds"))

ww_all = ww1 %>%
  # log
  dplyr::mutate(logvl=log(vl)) %>%
  mutate(vl=if_else(vl==0, 1, vl)) %>%
  # create indexes for INLA
  dplyr::mutate(day1=day,
                ara1=as.numeric(as.factor(ara_n)),
                ara2=ara1) %>% 
  # group KLBS and KLZH (otherwise not identifiable)
  dplyr::mutate(lab2=if_else(lab=="KLBS","KLZH",lab),
                lab_n2=as.numeric(as.factor(lab2))) %>% 
  # group lab and method
  dplyr::mutate(lab_method=factor(paste0(lab2,"_",method)),
                lab_method=relevel(lab_method, ref="EAWAG_0"),
                lab_method_n=as.numeric(lab_method),
                pop_totalb=pop_total/1000,
                prop_under_20b=u20*100,
                prop_over_65b=o65*100,
                prop_non_ch_eub=nec*100)
saveRDS(ww_all,file=paste0("../",controls$savepoint,"ww_all.rds"))

corr= readRDS(file=paste0("../",controls$savepoint,"corr_all_ara.rds"))

#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for deviation from average trend by ARA
# creation: jriou 
# init date: 2023-06-14
#:::::::::::::::::::::::::::::



library(patchwork)
library(tidyverse)
library(data.table)

model = readRDS(paste0("../",controls$savepoint,"ma5.3.2.rds"))

ww = readRDS(file=paste0("../",controls$savepoint,"ww_all.rds"))

corr= readRDS(file=paste0("../",controls$savepoint,"corr_all_ara.rds"))

#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for deviation from average trend by ARA
# creation: jriou 
# init date: 2023-06-14
#:::::::::::::::::::::::::::::

nara = length(unique(ww$ara1))
corr_periods = ww %>% 
  dplyr::select(date,period) %>% 
  dplyr::distinct()
corr_days = tibble(day=0:max(ww$day),
                   date=seq.Date(from=min(ww$date),to=max(ww$date),by=1)) %>% 
  left_join(corr_periods,by = join_by(date))
tt = model$summary.random$day1
alldays = unique(tt1$ID)
ndays = length(alldays)

tt = tt %>%  
  dplyr::bind_cols(day=rep(alldays,nara),
                   ara1=rep(1:nara,each=ndays)) %>% 
  dplyr::left_join(corr,by = join_by(ara1)) %>% 
  dplyr::left_join(corr_days,by="day") %>% 
  as_tibble()

  #tt2 = dplyr::filter(tt1,ara_name %in% extremes$ara_name) %>% 
  #  dplyr::mutate(ara_name=factor(ara_name,levels=extremes$ara_name)) 
  
  
tt_clusters = tibble()
  for(peri in 1:5){
    tt = tt1  %>% filter(period==peri)
    tt_masked = left_join(tt, ww %>% select(date, ara_name, vl), by=c('date', 'ara_name'))
    mask_in_period = tt_masked %>% 
      group_by(ara_id) %>%
      summarise(x = sum(vl, na.rm = TRUE))
    aras_to_include = mask_in_period %>% filter(x != 0) %>% select(ara_id)
    
    tt =  tt %>% filter((ara_id %in% unlist(aras_to_include)))
    
    g1 = tt %>% 
      ggplot() +
      geom_hline(yintercept=1,linetype=2,alpha=.5) +
      geom_ribbon(aes(x=date,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
      geom_line(aes(x=date,y=exp(mean),colour=NUTS2_name)) +
      facet_wrap(~ara_name) +
      scale_colour_discrete(guide="none") +
      scale_fill_discrete(guide="none") +
      scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
      coord_cartesian(ylim=c(.05,20)) +
      #labs(title=paste0("Top deviations from average time trend (top ",ntop,")"),x="Day",y="Relative viral load by ARA") + 
      theme(axis.text.x = element_text(angle=45,hjust=1))
    ggsave(g1, file=paste0('../outputs/plots/true_shapes_', peri, '.png'))
    
    tt_wide = data.table(dcast(tt, formula = day ~ ara_name, value.var = 'mean'))
    
    tt_list = as.list(tt_wide[,-c('day')])
    
    
    
    write.csv(tt, paste0('../data/true_shapes_' , peri, '.csv'))
    
    clustout_euclidean = dtwclust::tsclust(series=dtwclust::zscore(tt_list), 
                                           type='partitional', 
                                           distance = 'sbd',
                                           k=2:15, 
                                           return.objects=TRUE)
    
    names(clustout_euclidean) <- paste0(2:15)
    cvi_out = t(sapply(clustout_euclidean, dtwclust::cvi, type = "internal"))
    cvi_out = data.frame(cvi_out) %>% mutate(n_clusters = as.numeric(names(clustout_euclidean) ))
    cvi_out = cvi_out %>% pivot_longer(!n_clusters, names_to = 'type', values_to = 'score')
    cvi_out %>% 
      ggplot() +
      geom_line(aes(x=n_clusters, y=score, color=type))+
      facet_wrap(~type, scale='free')
    ggsave(file=paste0('../outputs/plots/cvi_euc_', peri, '.png'))
    
    
    
    clustout_tadpole = dtwclust::tsclust(series=dtwclust::zscore(tt_list), 
                                         type='tadpole',
                                         k=15, 
                                         control = dtwclust::tadpole_control(dc=3:15, window.size = 18L),
                                         return.objects=TRUE)
    
    names(clustout_tadpole) <- paste0(3:15)
    cvi_out = t(sapply(clustout_tadpole, dtwclust::cvi, type = "internal"))
    cvi_out = data.frame(cvi_out) %>% mutate(n_clusters = as.numeric(names(clustout_tadpole) ))
    cvi_out = cvi_out %>% pivot_longer(!n_clusters, names_to = 'type', values_to = 'score')
    cvi_out %>% 
      ggplot() +
      geom_line(aes(x=n_clusters, y=score, color=type))+
      facet_wrap(~type, scale='free')
    
    ggsave(file=paste0('../outputs/plots/cvi_tad_', peri, '_dc.png'))
    
    
    
    clustout_tadpole = dtwclust::tsclust(series=dtwclust::zscore(tt_list), 
                                         type='tadpole',
                                         k=2:25, 
                                         control = dtwclust::tadpole_control(dc=4, window.size = 18L),
                                         return.objects=TRUE)
    
    names(clustout_tadpole) <- paste0(2:25)
    cvi_out = t(sapply(clustout_tadpole, dtwclust::cvi, type = "internal"))
    cvi_out = data.frame(cvi_out) %>% mutate(n_clusters = as.numeric(names(clustout_tadpole) ))
    cvi_out = cvi_out %>% pivot_longer(!n_clusters, names_to = 'type', values_to = 'score')
    cvi_out %>% 
      ggplot() +
      geom_line(aes(x=n_clusters, y=score, color=type))+
      facet_wrap(~type, scale='free')
    
    ggsave(file=paste0('../outputs/plots/cvi_tad_', peri, '_k.png'))
    
    
    cluster_df = data.table(
      ara_name = names(tt_list)
    )
    
    cluster_df[, cluster_tad:=clustout_tadpole[[9]]@cluster]
    
    tt_cluster = merge(tt, cluster_df, by=c('ara_name'))
    
    
    g_clust = tt_cluster%>% 
      ggplot() +
      geom_hline(yintercept=1,linetype=2,alpha=.5) +
      #geom_ribbon(aes(x=date,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
      geom_line(aes(x=date,y=exp(mean), color=NUTS2_name, group=ara_name)) +
      facet_wrap(~cluster_tad) +
      scale_colour_discrete() +
      #scale_fill_discrete(guide="none") +
      scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
      coord_cartesian(ylim=c(.05,20)) +
      labs(title=paste0("Clusters"),x="Day",y="Relative viral load by ARA") + 
      theme(axis.text.x = element_text(angle=45,hjust=1))

    ggsave(g_clust, file=paste0('../outputs/plots/gclust_tad_', peri, '.png'))
    
    tt_cluster_sf = merge(shapes$ara_shp, unique(tt_cluster[,c('ara_id', 'ara_name', 'cluster_tad')], by=c('ara_id')))
    map = tt_cluster_sf %>% 
      ggplot() + 
      geom_sf(data=cantons) + 
      geom_sf(aes(fill = as.character(cluster_tad))) + 
      theme_minimal()
      
    ggsave(map, file=paste0('../outputs/plots/gclust_tad_map_', peri, '.png'))
    
    
    tt_cluster %>% mutate(period=peri)
    tt_clusters = rbind(tt_clusters, tt_cluster)
  }


### Transfer Entropy

tt_ind = 1:length(tt_list)

te = RTransferEntropy::transfer_entropy(x=unlist(tt_list[[6]]), y=unlist(tt_list[[7]]), quiet=T )

future::plan(future::multisession, workers = 6)
te_mat_list = purrr::map2(expand.grid(tt_ind, tt_ind)[,1], expand.grid(tt_ind, tt_ind)[,2], function(x,y){sum(RTransferEntropy::transfer_entropy(tt_list[[x]], tt_list[[y]], quiet=T)$coef[,1])})

kclustsim = cluster::pam(matrix(unlist(te_mat_list), nrow=length(tt_list)), 10, diss = T)

te_dt = rbind(tt_ind, unlist(te_mat_list))



matrix(unlist(te_mat_list), nrow=104)


shapes = readRDS(fs::path(controls$savepoint,"shapes.rds"))


catchments = shapes$ara_shp
catchments = sf::st_transform(catchments, 25830)
catchment_centroids = sf::st_centroid(catchments)

mtx_dist = sf::st_distance(catchment_centroids, catchment_centroids)

library(igraph)



kclustsim = cluster::pam(mtx_dist, 10, diss=TRUE)

spatial_cluster = kclustsim$clustering


spatial_clust = data.table::data.table(ara_id = catchment_centroids$ara_id, 
                              cluster_spat = spatial_cluster)


catchments_spat = merge(catchments, spatial_clust, by=c('ara_id'))


spa_clust_map_pm =
  ggplot() +  
  geom_sf(data=catchments_spat, aes(fill=as.character(cluster_spat)))


#stays = data.table::fread('../../02_data/other/stays.csv')


movements1  = data.table::fread('../../../02_data/other/movements_all_monthly/movements_2023-01-01--2023-01-31.csv')
movements2  = data.table::fread('../../../02_data/other/movements_all_monthly/movements_2023-02-01--2023-02-28.csv')
movements3  = data.table::fread('../../../02_data/other/movements_all_monthly/movements_2023-03-01--2023-03-31.csv')
movements4  = data.table::fread('../../../02_data/other/movements_all_monthly/movements_2023-04-01--2023-04-30.csv')
movements5  = data.table::fread('../../../02_data/other/movements_all_monthly/movements_2023-05-01--2023-05-31.csv')
movements6  = data.table::fread('../../../02_data/other/movements_all_monthly/movements_2023-06-01--2023-06-30.csv')
movements7  = data.table::fread('../../../02_data/other/movements_all_monthly/movements_2023-07-01--2023-07-31.csv')
movements8  = data.table::fread('../../../02_data/other/movements_all_monthly/movements_2023-08-01--2023-08-31.csv')
movements9  = data.table::fread('../../../02_data/other/movements_all_monthly/movements_2023-09-01--2023-09-30.csv')
movements10  = data.table::fread('../../../02_data/other/movements_all_monthly/movements_2023-10-01--2023-10-31.csv')
movements11  = data.table::fread('../../../02_data/other/movements_all_monthly/movements_2023-11-01--2023-11-30.csv')
movements12  = data.table::fread('../../../02_data/other/movements_all_monthly/movements_2023-12-01--2023-12-31.csv')

movements1[, uniques_01 := uniques]
movements2[, uniques_02 := uniques]
movements3[, uniques_03 := uniques]
movements4[, uniques_04 := uniques]
movements5[, uniques_05 := uniques]
movements6[, uniques_06 := uniques]
movements7[, uniques_07 := uniques]
movements8[, uniques_08 := uniques]
movements9[, uniques_09 := uniques]
movements10[, uniques_10 := uniques]
movements11[, uniques_11 := uniques]
movements12[, uniques_12 := uniques]

movements_long = rbind(
movements1, 
movements2, 
movements3, 
movements4, 
movements5, 
movements6, 
movements7, 
movements8, 
movements9, 
movements10,
movements11,
movements12)

movements_long[uniques=='<10', uniques:='10']
movements_long[is.na(uniques), uniques:=0]
movements_long[, uniques:=as.numeric(uniques)]

movements_wide = dcast(movements_long, formula = agent_type + origin_name + destination_name ~ month, value.var = 'uniques')


movements_mean = movements_long[, uniques_mean := mean(uniques, na.rm=T), by=c('agent_type', 'origin_name', 'destination_name')] 
movements_mean[, uniques:=uniques_mean]
movements_mean = movements_mean[, -c('uniques_mean', 'month')]

movements_mean = unique(movements_mean)

regular_movements = movements_mean[agent_type=='Regular' ,]
regular_movements = regular_movements[order(origin_name, destination_name),]


ww_all = data.table(ww_all)
ww_all[, ara_name_caps:=str_to_upper(ara_name)]
sort(unique(ww_all$ara_name))

id_names = unique(ww_all[, c('ara_id', 'ara_name_caps', 'NUTS2')])
id_names = id_names[!is.na(ara_name_caps), ]

regular_movements = merge(regular_movements, 
                          id_names, 
                          by.x = 'origin_name', 
                          by.y = 'ara_name_caps')

regular_movements[, origin_id:=ara_id]

regular_movements = regular_movements[, -c('ara_id')]

regular_movements = merge(regular_movements, 
                          id_names, 
                          by.x = 'destination_name', 
                          by.y = 'ara_name_caps')

regular_movements[, destination_id:=ara_id]

regular_movements = regular_movements %>% complete(origin_id, destination_id) %>% data.table()

regular_movements = regular_movements[order(origin_id, destination_id),]

regular_movements[, uniques_num := as.numeric(uniques)]



regular_movements[is.na(uniques_num), uniques_num:=1]
regular_movements[uniques_num < 1, uniques_num:=0]

mtx_mobil = matrix(regular_movements$uniques_num, nrow=length(unique(regular_movements$origin_id)))

(mtx_mobil + t(mtx_mobil))/2.0

mclustsim = cluster::pam( log(mtx_mobil), 10)

mobil_cluster = mclustsim$clustering


mobil_clusters = data.table::data.table(ara_id = sort(unique(regular_movements$origin_id)), 
                                        mobil_cluster = mobil_cluster)


catchments_mobil = merge(catchments, mobil_clusters, by=c('ara_id'))


mob_clust_map_pm =
  ggplot() +  
  geom_sf(data=catchments_mobil, aes(fill=as.character(mobil_cluster)))

mob_graph = graph_from_adjacency_matrix(mtx_mobil, mode = 'undirected', diag = F, weighted = T)

#mob_graph = graph_from_edgelist(as.matrix(regular_movements[,c('origin_name', 'destination_name')]))
#edge.attributes(mob_graph)$weight <- regular_movements$uniques_num

clusters = cluster_infomap(mob_graph)



mobil_cluster_ig = clusters$membership


mobil_clusters_ig = data.table::data.table(ara_id = sort(unique(regular_movements$origin_id)), 
                                        mobil_cluster = mobil_cluster_ig)


catchments_mobil_ig = merge(catchments, mobil_clusters_ig, by=c('ara_id'))



mob_clust_map_ig = 
  ggplot() +  
  geom_sf(data=catchments_mobil_ig, aes(fill=as.character(mobil_cluster)))


pal <- RColorBrewer::brewer.pal(length(unique(clusters$membership)), name ="Set3")

coords = sf::st_coordinates(catchment_centroids  %>% filter(ara_id %in% regular_movements$origin_id) %>% arrange(ara_id))
graph_plot_mob = 
  plot(mob_graph, 
     layout = as.matrix(coords), 
     edge.width = E(mob_graph)$weight/30000, 
     vertex.size=5, 
     vertex.label='', 
     edge.arrow.size= 0, 
     edge.curved = TRUE, 
     vertex.color=pal[clusters$membership])




tt_wide = data.table(dcast(tt, formula = day ~ ara_id, value.var = 'mean'))

tt_list = as.list(tt_wide[,-c('day')])

coords = sf::st_coordinates(catchment_centroids  %>% filter(ara_id %in% regular_movements$origin_id) %>% arrange(ara_id))

#coords = sf::st_coordinates(catchment_centroids  %>% filter(ara_id %in% tt$ara_id) %>% arrange(ara_id))
#dtw_graph = graph_from_adjacency_matrix(as.matrix(dtw_dists), mode = 'undirected', diag = F, weighted = T)
#clusters = cluster_infomap(dtw_graph)


nuts_clust_map = 
  ggplot() +  
  geom_sf(data=catchments_nuts, aes(fill=as.character(NUTS2)))



(nuts_clust_map + spa_clust_map_pm) / (mob_clust_map_pm + mob_clust_map_ig) +  plot_layout(guides = "collect")


igraph::centralization.closeness(mob_graph)

catchments_mobil_ig %>% mutate(eigen := eigen_centrality(mob_graph)$vector)

mob_graph = igraph::set.vertex.attribute(mob_graph, name = 'eigen', value = eigen_centrality(mob_graph)$vector) 
mob_graph = igraph::set.vertex.attribute(mob_graph, name = 'degree', value = strength(mob_graph)) 
mob_graph = igraph::set.vertex.attribute(mob_graph, name = 'page', value = page_rank(mob_graph)$vector)
mob_graph = igraph::set.vertex.attribute(mob_graph, name = 'between', value = betweenness(mob_graph, weights = 1.0/E(mob_graph)$weight)) 


features = data.table(st_drop_geometry(catchments_mobil_ig))

features = features %>% 
  mutate( 
    degree = V(mob_graph)$degree,
    betweenness = V(mob_graph)$between,
    e_centrality = V(mob_graph)$eigen,
    prank = V(mob_graph)$page
  )


high_deg = features %>% arrange(-degree) %>% head(12)

high_bet = features %>% arrange(-betweenness) %>% head(12)

high_cen = features %>% arrange(-e_centrality) %>% head(12)

high_pra = features %>% arrange(-prank) %>% head(12)



features = features %>% mutate(top_deg := ara_id %in% high_deg$ara_id)
features = features %>% mutate(top_bet := ara_id %in% high_bet$ara_id)
features = features %>% mutate(top_cen := ara_id %in% high_cen$ara_id)
features = features %>% mutate(top_pra := ara_id %in% high_pra$ara_id)



frame_colors <- c("lightgray", "limegreen")  # red for 0, blue for 1
frame_width <- c(1, 2)  # red for 0, blue for 1

jpeg("../outputs/plots/mob_select_graph.jpg", height=862, width=811)
graph_plots <- par(mfrow=c(2,2))
graph_plot_mob_deg = 
  plot(mob_graph, 
       vertex.frame.color = frame_colors[features$top_deg + 1],
       vertex.frame.width = frame_width[features$top_deg + 1],
       layout = as.matrix(coords), 
       edge.width = E(mob_graph)$weight/50000, 
       vertex.size=V(mob_graph)$degree/100000, 
       vertex.label='', 
       edge.arrow.size= 0, 
       edge.curved = TRUE, 
       vertex.color=pal[clusters$membership], 
       main='Degree')

graph_plot_mob_bet = 
  plot(mob_graph, 
       vertex.frame.color = frame_colors[features$top_bet + 1],
       vertex.frame.width = frame_width[features$top_bet + 1],
       vertex.frame.color='lightgray',
       layout = as.matrix(coords), 
       edge.width = E(mob_graph)$weight/50000, 
       vertex.size=V(mob_graph)$between/400, 
       vertex.label='', 
       edge.arrow.size= 0, 
       edge.curved = TRUE, 
       vertex.color=pal[clusters$membership], 
       main='Betweeness')

graph_plot_mob_eig = 
  plot(mob_graph, 
       vertex.frame.color = frame_colors[features$top_cen + 1],
       vertex.frame.width = frame_width[features$top_cen + 1],
       vertex.frame.color='lightgray',
       layout = as.matrix(coords), 
       edge.width = E(mob_graph)$weight/50000, 
       vertex.size=V(mob_graph)$eigen*20, 
       vertex.label='', 
       edge.arrow.size= 0, 
       edge.curved = TRUE, 
       vertex.color=pal[clusters$membership], 
       main='Eigenvalue')

graph_plot_mob_page = 
  plot(mob_graph, 
       vertex.frame.color = frame_colors[features$top_pra + 1],
       vertex.frame.width = frame_width[features$top_pra+ 1],
       vertex.frame.color='lightgray',
       layout = as.matrix(coords), 
       edge.width = E(mob_graph)$weight/50000, 
       vertex.size=V(mob_graph)$page*500, 
       vertex.label='', 
       edge.arrow.size= 0, 
       edge.curved = TRUE, 
       vertex.color=pal[clusters$membership], 
       main='Page rank')

dev.off()


write.csv(features, '../outputs/ara_selections_mobility.csv')





features_coms = data.table()

features[, idx := as.numeric(rownames(features))]

for(comm in 1:12){
  ids = features[mobil_cluster == comm,]$idx
  mob_graph_comm = igraph::induced_subgraph(mob_graph, ids)
  mob_graph_comm = igraph::set.vertex.attribute(mob_graph_comm, name = 'eigen', value = eigen_centrality(mob_graph_comm)$vector) 
  mob_graph_comm = igraph::set.vertex.attribute(mob_graph_comm, name = 'degree', value = strength(mob_graph_comm)) 
  mob_graph_comm = igraph::set.vertex.attribute(mob_graph_comm, name = 'page', value = page_rank(mob_graph_comm)$vector)
  mob_graph_comm = igraph::set.vertex.attribute(mob_graph_comm, name = 'between', value = betweenness(mob_graph_comm, weights = 1.0/E(mob_graph_comm)$weight)) 
  
  features_com = data.table(st_drop_geometry(catchments_mobil_ig))
  
  features_com = features_com[mobil_cluster == comm,] %>% 
    mutate( 
      degree = V(mob_graph_comm)$degree,
      betweenness = V(mob_graph_comm)$between,
      eigen = V(mob_graph_comm)$eigen,
      page = V(mob_graph_comm)$page
    )
  features_com[, mobil_cluster := comm]
  features_com[, top_deg := features_com$degree == max(features_com$degree) * 1.0]
  features_com[, top_bet := features_com$between == max(features_com$between)* 1.0]
  features_com[, top_cen := features_com$eigen == max(features_com$eigen)* 1.0]
  features_com[, top_pag := features_com$page == max(features_com$page)* 1.0]
  features_com[, idx := ids]
  features_coms = rbind(features_coms, features_com)

  
  }

features_coms = features_coms[order(as.numeric(idx)),]



frame_colors <- c("lightgray", "Purple")  # red for 0, blue for 1
frame_width <- c(1, 2)  # red for 0, blue for 1

jpeg("../outputs/plots/mob_select_graph_comms.jpg", height=862, width=811)
graph_plots <- par(mfrow=c(2,2))
graph_plot_mob_deg = 
  plot(mob_graph, 
       vertex.frame.color = frame_colors[features_coms$top_deg + 1],
       vertex.frame.width = frame_width[features_coms$top_deg + 1],
       layout = as.matrix(coords), 
       edge.width = E(mob_graph)$weight/50000, 
       vertex.size=features_coms$degree/50000, 
       vertex.label='', 
       edge.arrow.size= 0, 
       edge.curved = TRUE, 
       vertex.color=pal[features$mobil_cluster], 
       main='Degree')

graph_plot_mob_bet = 
  plot(mob_graph, 
       vertex.frame.color = frame_colors[features_coms$top_bet + 1],
       vertex.frame.width = frame_width[features_coms$top_bet + 1],
       vertex.frame.color='lightgray',
       layout = as.matrix(coords), 
       edge.width = E(mob_graph)$weight/50000, 
       vertex.size=features_coms$between/5, 
       vertex.label='', 
       edge.arrow.size= 0, 
       edge.curved = TRUE, 
       vertex.color=pal[clusters$membership], 
       main='Betweeness')

graph_plot_mob_eig = 
  plot(mob_graph, 
       vertex.frame.color = frame_colors[features_coms$top_cen + 1],
       vertex.frame.width = frame_width[features_coms$top_cen + 1],
       vertex.frame.color='lightgray',
       layout = as.matrix(coords), 
       edge.width = E(mob_graph)$weight/50000, 
       vertex.size=features_coms$eigen*10, 
       vertex.label='', 
       edge.arrow.size= 0, 
       edge.curved = TRUE, 
       vertex.color=pal[clusters$membership], 
       main='Eigenvalue')

graph_plot_mob_page = 
  plot(mob_graph, 
       vertex.frame.color = frame_colors[features_coms$top_pag + 1],
       vertex.frame.width = frame_width[features_coms$top_pag + 1],
       vertex.frame.color='lightgray',
       layout = as.matrix(coords), 
       edge.width = E(mob_graph)$weight/50000, 
       vertex.size=features_coms$page*30, 
       vertex.label='', 
       edge.arrow.size= 0, 
       edge.curved = TRUE, 
       vertex.color=pal[clusters$membership], 
       main='Page rank')

dev.off()


write.csv(features_coms, '../outputs/ara_selections_mobility_comms.csv')


features_coms[top_deg == T,]
features_coms[top_bet == T,]
features_coms[top_cen == T,]
features_coms[top_pag == T,]

catchments_clusters = catchments %>% merge(features_coms, by=c('ara_id'))

catchments_graph = catchments %>% merge(features, by=c('ara_id'))



deg_clust_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_clusters %>% filter(top_deg==T), aes(fill=as.character(mobil_cluster)), color=NA) +
  ggtitle('Degree cluster') + 
  scale_y_continuous(breaks = c()) + 
  theme(
    legend.position = 'right', 
    legend.title = element_text('Mobility communities')
  )

bet_clust_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_clusters %>% filter(top_bet==T), aes(fill=as.character(mobil_cluster))) +
  ggtitle('Betweeness cluster') + 
  theme(
    legend.position = 'none'
  )

cen_clust_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_clusters %>% filter(top_cen==T), aes(fill=as.character(mobil_cluster))) +
  ggtitle('Centrality cluster') + 
  theme(
    legend.position = 'none'
  )

pag_clust_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_clusters %>% filter(top_pag==T), aes(fill=as.character(mobil_cluster))) +
  ggtitle('Page rank cluster') + 
  theme(
    legend.position = 'none'
  )

deg_graph_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_graph %>% filter(top_deg==T), aes(fill=as.character(mobil_cluster))) +
  ggtitle('Degree graph') + 
  theme(
   legend.position = 'none'
  )

bet_graph_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_graph %>% filter(top_bet==T), aes(fill=as.character(mobil_cluster))) +
  ggtitle('Betweeness graph') + 
  theme(
    legend.position = 'none'
  )

cen_graph_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_graph %>% filter(top_cen==T), aes(fill=as.character(mobil_cluster))) +
  ggtitle('Centrality graph') + 
  theme(
    legend.position = 'none'
  )

pag_graph_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_graph %>% filter(top_pra==T), aes(fill=as.character(mobil_cluster))) +
  ggtitle('Page rank graph') + 
  theme(
    legend.position = 'none'
  )

layout = 
"
ABCD
EFGH"

deg_graph_map +  bet_graph_map +  cen_graph_map +  pag_graph_map + 
  deg_clust_map + bet_clust_map + cen_clust_map + pag_clust_map + 
  plot_layout(design = layout, guides = "collect") 


