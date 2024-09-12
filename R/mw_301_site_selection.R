library(patchwork)
library(tidyverse)
library(data.table)
library(dtwclust)
library(igraph)


catchments = shapes$ara_shp
catchments = sf::st_transform(catchments, 25830)
catchment_centroids = sf::st_centroid(catchments)

mtx_dist = sf::st_distance(catchment_centroids, catchment_centroids)



kclustsim = cluster::pam(mtx_dist, 10, diss=TRUE)

spatial_cluster = kclustsim$clustering


spatial_clust = data.table::data.table(ara_id = catchment_centroids$ara_id, 
                                       cluster_spat = spatial_cluster)


catchments_spat = merge(catchments, spatial_clust, by=c('ara_id'))


spa_clust_map_pm =
  ggplot() +  
  geom_sf(data=catchments_spat, aes(fill=as.character(cluster_spat)))


#stays = data.table::fread('../../02_data/other/stays.csv')


movements1  = data.table::fread('../../02_data/other/movements_all_monthly/movements_2023-01-01--2023-01-31.csv')
movements2  = data.table::fread('../../02_data/other/movements_all_monthly/movements_2023-02-01--2023-02-28.csv')
movements3  = data.table::fread('../../02_data/other/movements_all_monthly/movements_2023-03-01--2023-03-31.csv')
movements4  = data.table::fread('../../02_data/other/movements_all_monthly/movements_2023-04-01--2023-04-30.csv')
movements5  = data.table::fread('../../02_data/other/movements_all_monthly/movements_2023-05-01--2023-05-31.csv')
movements6  = data.table::fread('../../02_data/other/movements_all_monthly/movements_2023-06-01--2023-06-30.csv')
movements7  = data.table::fread('../../02_data/other/movements_all_monthly/movements_2023-07-01--2023-07-31.csv')
movements8  = data.table::fread('../../02_data/other/movements_all_monthly/movements_2023-08-01--2023-08-31.csv')
movements9  = data.table::fread('../../02_data/other/movements_all_monthly/movements_2023-09-01--2023-09-30.csv')
movements10  = data.table::fread('../../02_data/other/movements_all_monthly/movements_2023-10-01--2023-10-31.csv')
movements11  = data.table::fread('../../02_data/other/movements_all_monthly/movements_2023-11-01--2023-11-30.csv')
movements12  = data.table::fread('../../02_data/other/movements_all_monthly/movements_2023-12-01--2023-12-31.csv')

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


high_deg = features %>% arrange(-degree) %>% head(10)

high_bet = features %>% arrange(-betweenness) %>% head(10)

high_cen = features %>% arrange(-e_centrality) %>% head(10)

high_pra = features %>% arrange(-prank) %>% head(10)



features = features %>% mutate(top_deg := ara_id %in% high_deg$ara_id)
features = features %>% mutate(top_bet := ara_id %in% high_bet$ara_id)
features = features %>% mutate(top_cen := ara_id %in% high_cen$ara_id)
features = features %>% mutate(top_pra := ara_id %in% high_pra$ara_id)



frame_colors <- c("lightgray", "limegreen")  # red for 0, blue for 1
frame_width <- c(1, 2)  # red for 0, blue for 1

jpeg("outputs/plots/mob_select_graph.jpg", height=862, width=811)
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


write.csv(features, 'outputs/ara_selections_mobility.csv')





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

jpeg("outputs/plots/mob_select_graph_comms.jpg", height=862, width=811)
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


write.csv(features_coms, 'outputs/ara_selections_mobility_comms.csv')


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
  ggtitle('Commuting centers') + 
  scale_y_continuous(breaks = c()) + 
  theme(
    legend.position = 'right', 
    legend.title = element_text('Mobility communities')
  )

#bet_clust_map = 
#  ggplot() + 
#  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
#  geom_sf(data=catchments_clusters %>% filter(top_bet==T), aes(fill=as.character(mobil_cluster)), color=NA) +
#  ggtitle('Betweeness cluster') + 
#  theme(
#    legend.position = 'none'
#  )
#
#cen_clust_map = 
#  ggplot() + 
#  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
#  geom_sf(data=catchments_clusters %>% filter(top_cen==T), aes(fill=as.character(mobil_cluster)), color=NA) +
#  ggtitle('Centrality cluster') + 
#  theme(
#    legend.position = 'none'
#  )
#
#pag_clust_map = 
#  ggplot() + 
#  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
#  geom_sf(data=catchments_clusters %>% filter(top_pag==T), aes(fill=as.character(mobil_cluster)), color=NA) +
#  ggtitle('Page rank cluster') + 
#  theme(
#    legend.position = 'none'
#  )

deg_graph_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_graph %>% filter(top_deg==T), aes(fill=as.character(mobil_cluster)), color=NA) +
  ggtitle('Degree graph') + 
  theme(
    legend.position = 'none'
  )

bet_graph_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_graph %>% filter(top_bet==T), aes(fill=as.character(mobil_cluster)), color=NA) +
  ggtitle('Betweeness graph') + 
  theme(
    legend.position = 'none'
  )

cen_graph_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_graph %>% filter(top_cen==T), aes(fill=as.character(mobil_cluster)), color=NA) +
  ggtitle('Centrality graph') + 
  theme(
    legend.position = 'none'
  )

pag_graph_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_graph %>% filter(top_pra==T), aes(fill=as.character(mobil_cluster)), color=NA) +
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



layout = 
  "
AB
CD"

deg_graph_map +  bet_graph_map +  pag_graph_map + 
  deg_clust_map + 
  plot_layout(design = layout, guides = "collect") 

populations = unique(ww_all[,c('ara_id', 'pop_total')])[!is.na(pop_total),]

catchments_spat_dt = data.table(merge(catchments_spat, populations, by='ara_id'), all.x=T, all.y=F)

catchments_spat_dt[pop_total == -Inf, pop_total:= 0]
for(comm in 1:10){
  catchments_spat_dt[cluster_spat == comm, top_spat := catchments_spat_dt[cluster_spat == comm,]$pop_total == max(catchments_spat_dt[cluster_spat == comm,]$pop_total) * 1.0]
  
}
catchments_spat_dt[, top_spat := catchments_spat_dt$pop_total == max(catchments_spat_dt$pop_total) * 1.0]
spatial_ids = catchments_spat_dt[top_spat == TRUE,]$ara_id

spat_clust_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_spat %>% filter(ara_id %in% spatial_ids), aes(fill=as.character(cluster_spat)), color=NA) +
  ggtitle('Spatial clusters') + 
  theme(
    legend.position = 'none'
  )

nuts_frame = unique(ww_all[!is.na(NUTS2), c('ara_id', 'NUTS2', 'NUTS2_name', 'total_pop')])

for(comm in 1:7){
  nuts_frame[NUTS2 == comm, top_nuts := nuts_frame[NUTS2 == comm, ]$total_pop == max(nuts_frame[NUTS2 == comm, ]$total_pop)]
}

while(sum(nuts_frame$top_nuts) < 10){
  nuts_frame[top_nuts==F, top_nuts := nuts_frame[top_nuts==F, ]$total_pop == max(nuts_frame[top_nuts==F, ]$total_pop)]
}


nuts_ids = nuts_frame[top_nuts == T, ]$ara_id

nuts_pop_map = 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='gray')+
  geom_sf(data=catchments_spat %>% filter(ara_id %in% nuts_ids), aes(fill=as.character(cluster_spat)), color=NA) +
  ggtitle('Region based selection') + 
  theme(
    legend.position = 'none'
  )


layout = 
  "
AB
CD
EF"

nuts_pop_map + spat_clust_map + 
deg_graph_map +  bet_graph_map +  pag_graph_map + 
  deg_clust_map + 
  plot_layout(design = layout, guides = "collect") 




select_catchments = list(
  NUTS2 = nuts_ids, 
  spatial = spatial_ids, 
  mob_degree = (catchments_graph %>% filter(top_deg == T ))$ara_id, 
  mob_between = (catchments_graph %>% filter(top_bet == T ))$ara_id, 
  mob_pagerank = (catchments_graph %>% filter(top_pra == T))$ara_id, 
  mob_commute = unique((catchments_clusters %>% filter(top_deg == T))$ara_id)
  
)


saveRDS(object = select_catchments, file = 'outputs/select_catchments.rds')


select_shapes = shapes$ara_shp %>% mutate(selection = 0) %>% filter(ara_id == 0)
for(name in 1:length(selction_ids)){
  selection_n = shapes$ara_shp %>% filter(ara_id %in% selction_ids[[name]]) %>% mutate(selection=name)
  select_shapes = rbind(select_shapes, selection_n)
}


ggplot() + 
  geom_sf(data = shapes$canton_shp, alpha=0.2) + 
  geom_sf(data = shapes$ara_shp, alpha=0.2, color=NA, fill='lightseagreen') + 
  
  geom_sf(data = select_shapes, aes(fill=ara_id)) + 
  facet_wrap(~selection, labeller= labeller(selection = model_names))+
  scale_fill_discrete(guide=F)+
  theme_minimal() + 
  theme(panel.grid=element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())


