
source('setup.R')

## select model with daily deviations from GLOBAL trend 
model = readRDS(paste0("../",controls$savepoint,"ma5.4.2.rds"))

## load wastewater data
ww = readRDS(file=paste0("../",controls$savepoint,"ww_all.rds"))
corr = readRDS(file=paste0("../",controls$savepoint,"corr_all_ara.rds"))


## load shapes
shapes = readRDS(fs::path("../",controls$savepoint,"shapes.rds"))

## set number of aras with data 
nara = length(unique(ww$ara1))

## Generate 'deviation from global trend' time series
corr_periods = ww %>% 
  dplyr::select(date,period) %>% 
  dplyr::distinct()
corr_days = tibble(day=0:max(ww$day),
                   date=seq.Date(from=min(ww$date),to=max(ww$date),by=1)) %>% 
  left_join(corr_periods,by = join_by(date))
tt1 = model$summary.random$day1
alldays = unique(tt1$ID)
ndays = length(alldays)

tt1 = tt1 %>%  
  dplyr::bind_cols(day=rep(alldays,nara),
                   ara1=rep(1:nara,each=ndays)) %>% 
  dplyr::left_join(corr,by = join_by(ara1)) %>% 
  dplyr::left_join(corr_days,by="day") %>% 
  as_tibble()

  

## ------  cluster deviations for each period 1 - 5  ----------
## initialise containers
av_dists_all = data.table()
tt_cluster_all = data.table()
centroid_distances_all = data.table()


## iterate through periods
for (peri in 1:5){
  
  ## filter deviations for period selected
  tt = data.table(tt1  %>% filter(period %in% c(peri)))
  
  tt_masked = left_join(tt, ww %>% filter(period==peri) %>% select(date, ara_name, vl), by=c('date', 'ara_name'))
  
  ## remove aras with all NA values 
  mask_in_period = tt_masked %>% 
    group_by(ara_id) %>%
    summarise(x = sum(vl, na.rm = TRUE))
  
  mask_dates= tt_masked %>% 
    group_by(ara_id, day) %>%
    mutate(include = !is.na(vl))
  
  ## remove aras with more than 30% of samples missing
  mask_from_na = mask_dates %>% 
    group_by(ara_id) %>%
    mutate(x = sum(include, na.rm = TRUE)/length(include))
  
  aras_to_include_na = unique(mask_from_na %>% filter(x > 0.3) %>% select(ara_id))
  
  
  aras_to_include = mask_in_period %>% filter(x != 0) %>% select(ara_id)
  
  tt =  tt %>% filter((ara_id %in% unlist(aras_to_include)))
  tt =  tt %>% filter((ara_id %in% unlist(aras_to_include_na)))
  
  ## convert time-series to wide format and convert to list for clustering.
  tt_wide = data.table(dcast(tt, formula = day ~ ara_name, value.var = 'mean'))
  
  tt_list = as.list(tt_wide[,-c('day')])
  
  
  
  # run clustering algorithm for 2 to 10 clusters
  clustout_euclidean = dtwclust::tsclust(series=dtwclust::zscore(tt_list), 
                                         type='hierarchical', 
                                         distance = 'dtw2',
                                         k=2:10, 
                                         return.objects=TRUE)
  
<<<<<<< HEAD
  saveRDS(clustout_euclidean,paste0("../",controls$savepoint,"/dtw_outputs/cluster_", peri, '.rds'))
=======
  saveRDS(clustout_euclidean, paste0('../outputs/cluster_', peri, '.rds'))
>>>>>>> 363ac959dd96ba4efe9f05b2834a1646f4f35333
  
  ## evaluate partitions using CVI metrics - output as plot
  names(clustout_euclidean) <- paste0(2:10)
  cvi_out = t(sapply(clustout_euclidean, dtwclust::cvi, type = "valid"))
  cvi_out = data.frame(cvi_out) %>% mutate(n_clusters = as.numeric(names(clustout_euclidean) ))
  cvi_out = cvi_out %>% pivot_longer(!n_clusters, names_to = 'type', values_to = 'score')
  cvi_out %>% 
    ggplot() +
    geom_line(aes(x=n_clusters, y=score, color=type))+
    facet_wrap(~type, scale='free')
  
  ## Check max and mean cluster sizes at different n 
  plot(sapply(1:9, function(x){mean(clustout_euclidean[[x]]@clusinfo$size * clustout_euclidean[[x]]@clusinfo$av_dist/sum(clustout_euclidean[[x]]@clusinfo$size))})) 
  points(sapply(1:9, function(x){max(clustout_euclidean[[x]]@clusinfo$size * clustout_euclidean[[x]]@clusinfo$av_dist/sum(clustout_euclidean[[x]]@clusinfo$size))})) 
  
  
  for(n in 1:9){
    ## output the dtw distances betweeen deviations to estimate cluster strengths
<<<<<<< HEAD
    distmat = data.table::data.table(clustout_euclidean[[n]]@distmat)
    aras = rownames(clustout_euclidean[[n]]@distmat)
    distmat[, ara_1 := aras]
    distances = data.table::melt(distmat, id.vars = 'ara_1', variable.name = 'ara_2', value.name = 'Dist')
    
    clusters = data.table::data.table(ara = aras, cluster = clustout_euclidean[[n]]@cluster)
    
    distances_clust = merge(distances, clusters, by.x = 'ara_1', by.y = 'ara',by=.EACHI)
=======
    distmat = data.table(clustout_euclidean[[n]]@distmat)
    aras = rownames(clustout_euclidean[[n]]@distmat)
    distmat[, ara_1 := aras]
    distances = melt(distmat, id.vars = 'ara_1', variable.name = 'ara_2', value.name = 'Dist')
    
    clusters = data.table(ara = aras, cluster = clustout_euclidean[[n]]@cluster)
    
    distances_clust = merge(distances, clusters, by.x = 'ara_1', by.y = 'ara')
>>>>>>> 363ac959dd96ba4efe9f05b2834a1646f4f35333
    distances_clust = distances_clust[, cluster_1 := cluster][,-'cluster']
    distances_clust = merge(distances_clust, clusters, by.x = 'ara_2', by.y = 'ara')
    distances_clust = distances_clust[, cluster_2 := cluster][,-'cluster']
    
    ## Define strength of clusters using distances within and without clusters
    compare_distances = function(distances_clust, clust, ara){
      distances_to_common_clust = distances_clust[ara_1 == ara & ara_2 != ara & cluster_2 == clust, ]$Dist
      distances_to_other_clusts = distances_clust[ara_1 == ara  & cluster_2 != clust, ]$Dist
      quantile(distances_to_common_clust, probs = 0.5) < quantile(distances_to_other_clusts, probs = 0.05)
    }
    
    ## select strong clusters in data 
    clusters[, strong_clust := mapply(x=clusters$ara, y=clusters$cluster, FUN = function(x,y){compare_distances(distances_clust,y,x)}) ]
    
    ## merge clusters to deviations
    tt_cluster = merge(tt, clusters, by.x=c('ara_name'), by.y='ara')
    
    
    # identify the aras which represent the centroids of the clusters 
    cent_aras=names(clustout_euclidean[[n]]@centroids)
    # merged to spatial data and plot for records 
    tt_cluster_sf = merge(shapes$ara_shp, unique(tt_cluster[,c('ara_id', 'ara_name', 'cluster', 'strong_clust')], by=c('ara_id')))
    
    
    map= 
      ggplot() + 
      geom_sf(data=tt_cluster_sf, aes(fill = as.character(cluster)))+ 
      geom_sf(data=tt_cluster_sf %>% filter(ara_name %in% cent_aras), color='limegreen', linewidth=1, fill=NA)# + facet_wrap(~period)
    
    
    
<<<<<<< HEAD
    ggsave(paste0("../",controls$savepoint,'/dtw_plots/map_period_', peri, '_clust_',n, '.png'))
=======
    ggsave(paste0('../plots/map_period_', peri, '_clust_',n, '.png'))
>>>>>>> 363ac959dd96ba4efe9f05b2834a1646f4f35333
    
    
    map= 
      ggplot() + 
      geom_sf(data=tt_cluster_sf %>% filter(strong_clust==T), aes(fill = as.character(cluster))) + 
      geom_sf(data=tt_cluster_sf %>% filter(strong_clust==F), aes(fill = as.character(cluster)), alpha=0.2) + 
      geom_sf(data=tt_cluster_sf %>% filter(ara_name %in% cent_aras), color='limegreen', linewidth=1, fill=NA)# + facet_wrap(~period)
    
    
<<<<<<< HEAD
    ggsave(paste0("../",controls$savepoint,'/dtw_plots/map_period_strong_', peri, '_clust_',n, '.png'))
=======
    ggsave(paste0('../plots/map_period_strong_', peri, '_clust_',n, '.png'))
>>>>>>> 363ac959dd96ba4efe9f05b2834a1646f4f35333
    
    tt_cluster[, norm_mean := dtwclust::zscore(mean), by=c('ara_id')]
    
    
    ## plot time-series of clusters for records 
    g_clust = tt_cluster%>% #filter(strong_clust==T)%>% 
      ggplot() +
      geom_hline(yintercept=1,linetype=2,alpha=.5) +
      #geom_ribbon(aes(x=date,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
      geom_line(aes(x=date,y=exp(norm_mean), color=NUTS2_name, group=ara_name), alpha=0.5) +
      geom_line( data = tt_cluster%>%filter(ara_name %in% cent_aras), aes(x=date,y=exp(norm_mean), group=ara_name), color='limegreen', linetype=2, linewidth=2) +
      facet_wrap(~cluster) +
      scale_colour_discrete() +
      #scale_fill_discrete(guide="none") +
      scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
      coord_cartesian(ylim=c(.05,20)) +
      labs(title=paste0("Clusters"),x="Day",y="Relative viral load by ARA") + 
      theme(axis.text.x = element_text(angle=45,hjust=1))
    
    
    
<<<<<<< HEAD
    ggsave(paste0("../",controls$savepoint,'/dtw_plots/series_period_', peri, '_clust_',n, '.png'))
=======
    ggsave(paste0('../plots/series_period_', peri, '_clust_',n, '.png'))
>>>>>>> 363ac959dd96ba4efe9f05b2834a1646f4f35333
    
    g_clust = tt_cluster%>% filter(strong_clust==T)%>% 
      ggplot() +
      geom_hline(yintercept=1,linetype=2,alpha=.5) +
      #geom_ribbon(aes(x=date,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
      geom_line(aes(x=date,y=exp(norm_mean), color=NUTS2_name, group=ara_name), alpha=0.5) +
      geom_line( data = tt_cluster%>%filter(ara_name %in% cent_aras), aes(x=date,y=exp(norm_mean), group=ara_name), color='limegreen', linetype=2, linewidth=2) +
      facet_wrap(~cluster) +
      scale_colour_discrete() +
      #scale_fill_discrete(guide="none") +
      scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
      coord_cartesian(ylim=c(.05,20)) +
      labs(title=paste0("Clusters"),x="Day",y="Relative viral load by ARA") + 
      theme(axis.text.x = element_text(angle=45,hjust=1))
    
    
    
<<<<<<< HEAD
    ggsave(paste0("../",controls$savepoint,'/dtw_plots/series_period_strong_', peri, '_clust_',n, '.png'))
=======
    ggsave(paste0('../plots/series_period_strong_', peri, '_clust_',n, '.png'))
>>>>>>> 363ac959dd96ba4efe9f05b2834a1646f4f35333
    
    
    av_dists = data.table(clustout_euclidean[[n]]@clusinfo) 
    colnames(av_dists) = c('No. in cluster', 'Average distance')
    av_dists[, cluster := 1:nrow(av_dists)]
    av_dists[, period := peri]
    av_dists[, nclust := ..n]
    av_dists_all = rbind(av_dists_all, av_dists)
    
    tt_cluster[, period:=peri]
    tt_cluster[, nclust:=..n]
    tt_cluster_all = rbind(tt_cluster_all, tt_cluster)
    centroid_distances = data.table()
    for(ara1 in cent_aras){
      for(ara2 in cent_aras){
        centroid_distances_pair = data.table(
          ara_1 = c(ara1), 
          ara_2 = c(ara2), 
          cluster_1 = clusters[ara==ara1]$cluster,
          cluster_2 = clusters[ara==ara2]$cluster,
          distance = c(dtwclust::dtw2(tt_list[ara1], tt_list[ara2]))$distance, 
          period = c(peri), 
          nclust = n)
        centroid_distances = rbind(centroid_distances, centroid_distances_pair)
        
      }
      centroid_distances_all = rbind(centroid_distances_all, centroid_distances)
    }
  }
  
}

tt_cluster_all
av_dists_all %>% ggplot() + 
  ggbeeswarm::geom_beeswarm(aes(x =period, y=`Average distance`, color=period)) + 
  facet_wrap(~nclust)

centroid_distances_all[!(ara_1 == ara_2)] %>% ggplot() + 
  geom_point(aes(x = rank(distance), y=`distance`, color=as.character(period))) + 
  facet_wrap(~nclust)


centroid_distances_all_av = merge(centroid_distances_all, av_dists_all[,c('Average distance', 'cluster', 'period', 'nclust')], by.x = c('cluster_1', 'period', 'nclust'), by.y = c('cluster', 'period', 'nclust'))


centroid_distances_all_av[, rat_av_cen_dist := distance/`Average distance`]


centroid_distances_all_av[!(ara_1 == ara_2)] %>% ggplot() + 
  geom_point(aes(x = rank(rat_av_cen_dist), y=rat_av_cen_dist, color=as.character(period))) + 
  facet_wrap(~nclust)


centroid_distances_all_av[!(ara_1 == ara_2) & !is.infinite(rat_av_cen_dist)] %>% ggplot() + 
  geom_violin(aes(x = period, y=rat_av_cen_dist, color=as.character(period))) + 
  geom_point(aes(x = period, y=rat_av_cen_dist) ) + 
  facet_wrap(~nclust)


write.csv(centroid_distances_all_av, '../outputs/dtw_results/centroid_distances_all_av.csv')

write.csv(av_dists_all, '../outputs/dtw_results/av_dists_all.csv')

write.csv(tt_cluster_all, '../outputs/dtw_results/tt_cluster_all.csv')


# Maps of all nclust at all phases 

tt_cluster_sf_all = merge(shapes$ara_shp, unique(tt_cluster_all[,c('ara_id', 'ara_name', 'cluster', 'strong_clust', 'period', 'nclust')]), by=c('ara_id'))
map_all= 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='darkgray', linewidth=0.5, alpha=0.5) +
  geom_sf(data=tt_cluster_sf_all, aes(fill = as.character(cluster)), color=NA, alpha=0.8) + 
  facet_grid(nclust~period) + 
  geom_sf(data=shapes$see_shp, fill='midnightblue', color=NA)+ 
  theme_minimal() + 
  scale_fill_manual(name='Cluster',  values = paletteer_d("ggthemes::hc_darkunica")) +
  ggtitle('A') +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank())


ggsave(paste0('../plots/map_all.png'), plot = map_all)

map_all_strong= 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='darkgray', linewidth=0.5, alpha=0.5) +
  geom_sf(data=tt_cluster_sf_all %>% filter(strong_clust==T), aes(fill = as.character(cluster)), color=NA, alpha=0.8) + 
  geom_sf(data=tt_cluster_sf_all %>% filter(strong_clust==F), aes(color = as.character(cluster)), fill=NA, linewidth=0.1) + 
  facet_grid(nclust~period) + 
  geom_sf(data=shapes$see_shp, fill='midnightblue', color=NA)+ 
  theme_minimal() + 
  scale_fill_manual(name='Cluster',  values = paletteer_d("ggthemes::hc_darkunica")) +
  scale_color_discrete(guide=F)+
  scale_x_continuous(breaks = c())+
  ggtitle('B')+
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank()
  )


ggsave(paste0('../plots/map_all_strong_n.png'), plot=map_all_strong, height = 15, width=12)



map_all_choice= 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='darkgray', linewidth=0.5, alpha=0.5) +
  geom_sf(data=tt_cluster_sf_all %>% filter(nclust==7), aes(fill = as.character(cluster)), color=NA, alpha=0.8) + 
  facet_grid(nclust~period) + 
  geom_sf(data=shapes$see_shp, fill='midnightblue', color=NA)+ 
  theme_minimal() + 
  scale_fill_discrete(name='Cluster') + 
  ggtitle('A') +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank())


ggsave(paste0('../plots/map_all_choice.png'), plot = map_all_choice)

map_all_strong_choice= 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='darkgray', linewidth=0.5, alpha=0.5) +
  geom_sf(data=tt_cluster_sf_all %>% filter(nclust==7, strong_clust==T), aes(fill = as.character(cluster)), color=NA, alpha=0.8) + 
  geom_sf(data=tt_cluster_sf_all %>% filter(nclust==7, strong_clust==F), aes(color = as.character(cluster)), fill=NA, linewidth=0.1) + 
  facet_grid(nclust~period) + 
  geom_sf(data=shapes$see_shp, fill='midnightblue', color=NA)+ 
  theme_minimal() + 
  scale_fill_discrete(name='Cluster') + 
  scale_color_discrete(guide=F)+
  scale_x_continuous(breaks = c())+
  ggtitle('B')+
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank()
  )


ggsave(paste0('../plots/map_all_strong_choice.png'), plot=map_all_strong_choice, height = 15, width=12)


#pal <- RColorBrewer::brewer.pal(8, name ="Set3")

#for(peri in 1:5) {
#  centroid_graph = igraph::graph_from_data_frame(unique(centroid_distances_all[period==peri & distance > 0,c('cluster_1', 'cluster_2', 'distance')]), directed = F)
#  igraph::edge.attributes(centroid_graph)$weight <- log(1.0/(unique(centroid_distances_all[period==peri & distance > 0,c('cluster_1', 'cluster_2', 'distance')])$distance))
#  igraph::E(centroid_graph)$weight
#  plot(centroid_graph, 
#       edge.arrow.size= 0, 
#       edge.curved = TRUE, 
#       edge.width = igraph::E(centroid_graph)$weight * 10,
#       vertex.color=pal[ as.numeric(igraph::V(centroid_graph)$name)],
#       main=paste0('period ' ,peri))
#}

#image(as.matrix(dcast(unique(centroid_distances_all[period==peri, c('cluster_1', 'cluster_2', 'distance')]), formula=cluster_1~cluster_2, value.var = 'distance')))


# charts of select clusters based on maintaining geographic structure 

# relabling clusters to maintain similarity between phases 
tt_cluster_select = data.table(tt_cluster_all %>% filter((period==1 & nclust==7) | (period==2 & nclust==3) | (period==3 & nclust==2) | (period==4 & nclust == 2) | (period == 5 & nclust==3) ) )

cluster_swaps = find_cluster_swaps(tt_cluster_select)
tt_cluster_select[period==1, reassigned_cluster := cluster]
for(peri in 2:5){
  for(clust in 1:length(cluster_swaps[[1]][[peri]]))
    tt_cluster_select[period==peri & cluster==cluster_swaps[[1]][[peri]][[clust]], reassigned_cluster := cluster_swaps[[2]][[peri]][[clust]]]
}

# distances between centroids of the clusters
unique_distances = unique(centroid_distances_all[!(cluster_1 == cluster_2) , c('cluster_1', 'cluster_2', 'distance', 'period', 'nclust')]) 


unique_distances[cluster_1 < cluster_2, distance:=NA]

cluster_distance_mat = unique_distances[!is.na(distance) & nclust==7,] %>% ggplot() + 
  geom_tile(aes(x=cluster_1, y=cluster_2, fill=distance), )+
  geom_text(aes(x=cluster_1, y=cluster_2,label = round(distance,1)), color = "white", size = 3.7) +
  scale_fill_viridis_c(trans='reverse', guide=F)+
  scale_y_continuous(breaks=1:8)+
  scale_x_continuous(breaks=1:8)+
  facet_wrap(~period, nrow=1) + 
  ylab('Cluster')+
  xlab('Cluster')+
  ggtitle('C') + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
  )


cluster_distance_mat_1 = unique_distances[!is.na(distance) & ((period==1 & nclust==7)),] %>% ggplot() + 
  geom_tile(aes(x=cluster_1, y=cluster_2, fill=distance), )+
  geom_text(aes(x=cluster_1, y=cluster_2,label = round(distance,1)), color = "white", size = 3.7) +
  scale_fill_viridis_c(guide=F, limits=c(25,0), trans='reverse')+
  scale_y_continuous(breaks=1:8)+
  scale_x_continuous(breaks=1:8)+
  ylab('Cluster')+
  xlab('Cluster')+
  ggtitle('Phase 1') + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
  )

cluster_distance_mat_2 = unique_distances[!is.na(distance) & ( (period==2 & nclust== 3)),] %>% ggplot() + 
  geom_tile(aes(x=cluster_1, y=cluster_2, fill=distance), )+
  geom_text(aes(x=cluster_1, y=cluster_2,label = round(distance,1)), color = "white", size = 3.7) +
  scale_fill_viridis_c(guide=F, limits=c(25,0), trans='reverse')+
  scale_y_continuous(breaks=1:4, labels=cluster_swaps[[2]][[2]][order(cluster_swaps[[1]][[2]])])+
  scale_x_continuous(breaks=1:4, labels=cluster_swaps[[2]][[2]][order(cluster_swaps[[1]][[2]])])+
  ylab('Cluster')+
  xlab('Cluster')+
  ggtitle('Phase 2') + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
  )


cluster_distance_mat_3 = unique_distances[!is.na(distance) & ((period==3 & nclust==2) ),] %>% ggplot() + 
  geom_tile(aes(x=cluster_1, y=cluster_2, fill=distance), )+
  geom_text(aes(x=cluster_1, y=cluster_2,label = round(distance,1)), color = "white", size = 3.7) +
  scale_fill_viridis_c(guide=F, limits=c(25,0), trans='reverse')+
  scale_y_continuous(breaks=1:3, labels=cluster_swaps[[2]][[3]][order(cluster_swaps[[1]][[3]])])+
  scale_x_continuous(breaks=1:3, labels=cluster_swaps[[2]][[3]][order(cluster_swaps[[1]][[3]])])+
  facet_wrap(~period, nrow=1, scale='free') + 
  ylab('Cluster')+
  xlab('Cluster')+
  ggtitle('Phase 3') + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
  )


cluster_distance_mat_4 = unique_distances[!is.na(distance) & ( (period==4 & nclust==2)),] %>% ggplot() + 
  geom_tile(aes(x=cluster_1, y=cluster_2, fill=distance), )+
  geom_text(aes(x=cluster_1, y=cluster_2,label = round(distance,1)), color = "white", size = 3.7) +
  scale_fill_viridis_c(guide=F, limits=c(25,0), trans='reverse')+
  scale_y_continuous(breaks=1:3, labels=cluster_swaps[[2]][[4]][order(cluster_swaps[[1]][[4]])])+
  scale_x_continuous(breaks=1:3, labels=cluster_swaps[[2]][[4]][order(cluster_swaps[[1]][[4]])])+
  facet_wrap(~period, nrow=1, scale='free') + 
  ylab('Cluster')+
  xlab('Cluster')+
  ggtitle('Phase 4') + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
  )

cluster_distance_mat_5 = unique_distances[!is.na(distance) & ((period==5 & nclust==3)),] %>% ggplot() + 
  geom_tile(aes(x=cluster_1, y=cluster_2, fill=distance), )+
  geom_text(aes(x=cluster_1, y=cluster_2,label = round(distance,1)), color = "white", size = 3.7) +
  scale_fill_viridis_c(guide=F, limits=c(25,0), trans='reverse')+
  scale_y_continuous(breaks=1:4, labels=cluster_swaps[[2]][[5]][order(cluster_swaps[[1]][[5]])])+
  scale_x_continuous(breaks=1:4, labels=cluster_swaps[[2]][[5]][order(cluster_swaps[[1]][[5]])])+
  facet_wrap(~period, nrow=1, scale='free') + 
  ylab('Cluster')+
  xlab('Cluster')+
  ggtitle('Phase 5') + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
  )

centroid_dist_plot = cluster_distance_mat_1 + cluster_distance_mat_2 + cluster_distance_mat_3 + cluster_distance_mat_4 + cluster_distance_mat_5 + plot_layout(design = 'abcde')

ggsave(filename = '../plots/centroid_dist_plot.pdf', centroid_dist_plot, width=15, height=4)

joint_plot = map_all/cluster_distance_mat + plot_layout(heights=c(5,2))

ggsave(paste0('../plots/map_all_distances.png'), width=20, height=10, plot = joint_plot)

joint_plot_sc = map_all_strong/cluster_distance_mat + plot_layout(heights=c(5,2))

ggsave(paste0('../plots/map_all_distances_sc.png'), width=20, height=10, plot = joint_plot_sc)

joint_plot_all = map_all/map_all_strong/cluster_distance_mat + plot_layout(heights=c(4,4,3), guides = "collect") & theme(legend.position = 'bottom')
ggsave(paste0('../plots/map_all_distances_all_m.png'), width=12, height=8, plot = joint_plot_all)

ggsave(paste0('../plots/map_all_distances_all_n.pdf'), width=12, height=8, plot = joint_plot_all)



g_clust = tt_cluster_all %>% filter(nclust %in% 2:7 ) %>% 
  ggplot() +
  geom_hline(yintercept=1,linetype=2,alpha=.5) +
  #geom_ribbon(aes(x=date,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
  geom_line(aes(x=date,y=exp(norm_mean), color=NUTS2_name, group=ara_name), alpha=0.5) +
  #geom_line( data = tt_cluster_all %>% filter(period==1 & ara_name %in% cent_aras), aes(x=date,y=exp(norm_mean), group=ara_name), color='limegreen', linetype=2, linewidth=2) +
  facet_grid(cluster ~ nclust + period, scales = 'free') +
  scale_colour_discrete() +
  #scale_fill_discrete(guide="none") +
  scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
  coord_cartesian(ylim=c(.05,20)) +
  labs(title=paste0("Clusters"),x="Day",y="Relative viral load by ARA") + 
  theme(axis.text.x = element_text(angle=90,hjust=1))

png(filename = '../plots/all_series.png', width = 40, heigh=12, units='cm', res = 300) # open the file
g_clust
dev.off() # close the file


g_clust_png <- image_read('../plots/all_series.png')

# rotate
image_rotate(g_clust_png, 90) %>% image_write('../plots/all_series_rot.png')




g_clust = tt_cluster_all %>% filter(nclust==8 & strong_clust) %>% 
  ggplot() +
  geom_hline(yintercept=1,linetype=2,alpha=.5) +
  #geom_ribbon(aes(x=date,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
  geom_line(aes(x=date,y=exp(norm_mean), color=NUTS2_name, group=ara_name), alpha=0.5) +
  #geom_line( data = tt_cluster_all %>% filter(period==1 & ara_name %in% cent_aras), aes(x=date,y=exp(norm_mean), group=ara_name), color='limegreen', linetype=2, linewidth=2) +
  facet_grid(cluster~period, scales = 'free') +
  scale_colour_discrete() +
  #scale_fill_discrete(guide="none") +
  scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
  coord_cartesian(ylim=c(.05,20)) +
  labs(title=paste0("Clusters"),x="Day",y="Relative viral load by ARA") + 
  theme(axis.text.x = element_text(angle=45,hjust=1))


for (peri in 1:5){
  clust = readRDS(paste0('../outputs/cluster_', peri, '.rds'))
  dend = as.dendrogram(clust[[9]])

  nuts = unique(data.table(tt1)[ara_name %in% clust[[9]]$labels, c('ara_name', 'NUTS2_name')][order(ara_name),])$NUTS2_name
  colors_to_use <- as.numeric(as.factor(nuts))
  colors_to_use
  # But sort them based on their order in dend:
  colors_to_use <- colors_to_use[order.dendrogram(dend)]
  colors_to_use
  
  pdf(paste0('../plots/dendogram_',peri,'.pdf'), width = 15, height=17 )
  # Now we can use them
  dendextend::labels_colors(dend) <- gg_color_hue(7)[colors_to_use]
  dendextend::labels_cex(dend) = 1
  par(mar = c(5,1,1,5))
  plot(dend, horiz = TRUE) 
  
 
  dev.off()
  
}




tt_cluster_select[, cluster:=reassigned_cluster]
tt_cluster_sf_select = merge(shapes$ara_shp, unique(tt_cluster_select[,c('ara_id', 'ara_name', 'cluster', 'strong_clust', 'period', 'nclust')]), by=c('ara_id'))


av_cluster_strength[period==1, reassigned_cluster := cluster]
for(peri in 2:5){
  for(clust in 1:length(cluster_swaps[[1]][[peri]]))
    av_cluster_strength[period==peri & cluster==cluster_swaps[[1]][[peri]][[clust]], reassigned_cluster := cluster_swaps[[2]][[peri]][[clust]]]
}

av_cluster_strength = av_cluster_strength[, Cluster:=as.integer(reassigned_cluster)][, -c('reassigned_cluster')]
av_cluster_strength = av_cluster_strength[, Phase:=period][, -c('reassigned_cluster')]
av_cluster_strength[, `Normalised average distance` := `Average distance` /(`No. in cluster`**0.5)]

print(xtable::xtable(av_cluster_strength[, c('Phase', 'Cluster', 'No. in cluster', 'Average distance', 'Normalised average distance' )], 
                     caption = "Size of and average distance between member WWTPs of clusters in each phases of the analyisis"),
      include.rownames = F, 
      hline.after=c(-1, 0, cumsum(unique(av_cluster_strength[, c('period', 'nclust')])$nclust+1)), 
      caption.placement = "top"
      )
      


names_periods = paste0('Period ', 1:5 )
names(names_periods) = 1:5

mixed_map = tt_cluster_sf_select %>% filter(period %in% c(1,2,5) ) %>% ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='darkgray', linewidth=0.5, alpha=0.8) +
  geom_sf(aes(fill = as.character(cluster)), color=NA, alpha=0.8) + 
  scale_fill_manual(name='Cluster',  values = paletteer_d("ggthemes::hc_darkunica")) +
  geom_sf(data=shapes$see_shp, fill='midnightblue', color=NA)+ 
  theme_minimal() + 
  facet_wrap(~period, nrow=1, labeller = labeller(period = names_periods)) + 
  scale_color_discrete(guide=F)+
  scale_x_continuous(breaks = c())+
  ggtitle('A')+
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank()
  )
  
ggsave(paste0('../plots/map_mixed_clustsize.png'), width=20, height=10, plot = mixed_map)

names_clusters = paste0('Cluster ', 1:10 )
names(names_clusters) = 1:10


gclust_1= tt_cluster_select %>% filter((period==1 & nclust==7) ) %>% #| (period==2 & nclust==5) | (period == 5 & nclust==3) ) %>% 
  ggplot() +
  geom_hline(yintercept=1,linetype=2,alpha=.5) +
  #geom_ribbon(aes(x=date,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
  geom_line(aes(x=as.Date(date),y=exp(norm_mean), color=NUTS2_name, group=ara_name), alpha=0.8) +
  #geom_line( data = tt_cluster_all %>% filter(period==1 & ara_name %in% cent_aras), aes(x=date,y=exp(norm_mean), group=ara_name), color='limegreen', linetype=2, linewidth=2) +
  facet_wrap(~cluster, ncol=2, labeller = labeller(cluster = names_clusters)) +
  scale_color_discrete(name='Region') +    #scale_fill_discrete(guide="none") +
  scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
  scale_x_date(date_breaks='1 month', date_labels = '%b %y')+
  coord_cartesian(ylim=c(.05,20)) +
  labs(title=paste0("B"),x="",y="Relative viral load by ARA") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45,hjust=1, size=8)) 
  


gclust_2 = tt_cluster_select %>% filter((period==2 & nclust==3)) %>% #| (period==2 & nclust==5) | (period == 5 & nclust==3) ) %>% 
  ggplot() +
  geom_hline(yintercept=1,linetype=2,alpha=.5) +
  #geom_ribbon(aes(x=date,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
  geom_line(aes(x=as.Date(date),y=exp(norm_mean), color=NUTS2_name, group=ara_name), alpha=0.8) +
  #geom_line( data = tt_cluster_all %>% filter(period==1 & ara_name %in% cent_aras), aes(x=date,y=exp(norm_mean), group=ara_name), color='limegreen', linetype=2, linewidth=2) +
  facet_wrap(~cluster, ncol=2, labeller = labeller(cluster = names_clusters)) +
  scale_color_discrete(name='Region') +  
  scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
  scale_x_date(date_breaks='1 month', date_labels = '%b %y')+
  coord_cartesian(ylim=c(.05,20)) +
  labs(title=paste0("C"),x="",y="") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45,hjust=1, size=8)) 

gclust_5 = tt_cluster_select %>% filter((period == 5 & nclust==3)) %>% #| (period==2 & nclust==5) | (period == 5 & nclust==3) ) %>% 
  ggplot() +
  geom_hline(yintercept=1,linetype=2,alpha=.5) +
  #geom_ribbon(aes(x=date,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
  geom_line(aes(x=as.Date(date),y=exp(norm_mean), color=NUTS2_name, group=ara_name), alpha=0.8) +
  #geom_line( data = tt_cluster_all %>% filter(period==1 & ara_name %in% cent_aras), aes(x=date,y=exp(norm_mean), group=ara_name), color='limegreen', linetype=2, linewidth=2) +
  facet_wrap(~cluster, ncol=2, labeller = labeller(cluster = names_clusters)) +
  scale_color_discrete(name='Region') +  #scale_fill_discrete(guide="none") +
  scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
  scale_x_date(date_breaks='1 month', date_labels = '%b %y')+
  coord_cartesian(ylim=c(.05,20)) +
  labs(title=paste0("D"),x="",y="") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45,hjust=1, size=8)) 



saveRDS(object = tt_cluster_sf_all, file = '../outputs/tt_cluster_sf_all.rds')

library(ggsankeyfier)


clust_select_wide = dcast(unique(tt_cluster_select[, c('period', 'ara_name', 'cluster')]), 
                          formula =  ara_name ~ period, 
                          value.var = 'cluster')


colnames(clust_select_wide) = c('ara_name', paste0('Phase ', 1:5))

flow_df = clust_select_wide %>% make_long(`Phase 1`,  `Phase 2`,  `Phase 3`,  `Phase 4`, `Phase 5`)

sankeyflow = ggplot(flow_df %>% filter(!is.na(node)),  aes(x = x, 
                    next_x = next_x, 
                    node = node, 
                    next_node = next_node,
                    fill = factor(node),
                    label = node), alpha=0.5) +
  geom_sankey(na.rm=T, alpha=0.5, color=NA, width = 0.0005, space = 40) +
  geom_sankey_label(na.rm = T, width = 0.1, color=NA, space = 40, alpha=0.8) +
  geom_sankey_text(color='white', space = 40)+
  ggtitle('E') +
  scale_fill_manual(name='Cluster',  values = paletteer_d("ggthemes::hc_darkunica")) +
  theme_sankey(base_size = 6) + 
  xlab('') + 
  theme(legend.position = 'none', 
        plot.margin=margin(0,0,0,0),
        axis.text.x = element_text(angle=45,hjust=1, size=8), 
        title = element_text(size=12))


clust_select_wide




library(ggalluvial)

ggplot(data = clust_select_wide,
       aes(axis1 = `Phase 1`, axis2 = `Phase 2`, axis3=`Phase 3`, axis4=`Phase 4`, axis5=`Phase 5`)) +
  geom_alluvium(aes(fill=value)) +
  geom_stratum(na.rm=T) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void() + 
  scale_fill_viridis_c()

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cluster_flow_plot = 
  ggplot(unique(tt_cluster_select[, c('period', 'ara_name', 'cluster')]),
       aes(x = period, y=ara_name , stratum = as.character(cluster), alluvium = ara_name,
           fill = as.character(cluster), label = as.character(cluster))) +
  scale_x_discrete(expand = c(.1, .1)) + 
  
  geom_flow() +
  geom_stratum(color = "white", width = 1/6, linewidth=1) +
  geom_text(stat = "stratum", size = 2) +
  scale_fill_manual(name='Cluster',  values = paletteer_d("ggthemes::hc_darkunica")) +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 8, face = "bold"), 
    legend.position = 'none', 
    plot.margin=margin(0,0,0,0)
  )


ggsave(filename = '../plots/cluster_flow_plot.png', plot = cluster_flow_plot, width=12, height=4)
ggsave(filename = '../plots/cluster_flow_plot.pdf', plot = cluster_flow_plot, width=12, height=4)

cluster_flow_plot_strong = 
  ggplot(unique(tt_cluster_select[strong_clust==T, c('period', 'ara_name', 'cluster')]),
         aes(x = period, y=ara_name , stratum = as.character(cluster), alluvium = ara_name,
             fill = as.character(cluster), label = as.character(cluster))) +
  scale_x_discrete(expand = c(.1, .1)) + 
  scale_fill_manual(name='Cluster',  values = paletteer_d("ggthemes::hc_darkunica")) +
  geom_flow() +
  geom_stratum(color = "white", width = 1/6, linewidth=1) +
  geom_text(stat = "stratum", size = 4) +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", 
    )
  )


ggsave(filename = '../plots/cluster_flow_plot_strong.png', plot = cluster_flow_plot_strong, width=12, height=8)
ggsave(filename = '../plots/cluster_flow_plot_strong.pdf', plot = cluster_flow_plot_strong, width=12, height=8)


layout=
  "
aaaaaa
aaaaaa
aaaaaa
bbccdd
bbccdd
bbccdd
bbeeee
bbeeee
bbeeee
"


mixed_map_series = mixed_map + gclust_1 + gclust_2 + gclust_5 + sankeyflow + plot_layout(design=layout, guides = 'collect')

ggsave(paste0('../plots/map_series_mixed_clustsize.png'), width=10, height=9, plot = mixed_map_series)
ggsave(paste0('../plots/map_series_mixed_clustsize.pdf'), width=10, height=9, plot = mixed_map_series)

