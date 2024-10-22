library(patchwork)
library(tidyverse)
library(data.table)

## select model with daily deviations from GLOBAL trend 
model = readRDS(paste0("../",controls$savepoint,"ma5.4.4.rds"))

## load wastewater data
ww = readRDS(file=paste0("../",controls$savepoint,"ww_all.rds"))
corr = readRDS(file=paste0("../",controls$savepoint,"corr_all_ara.rds"))

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
  
  ## evaluate partitions using CVI metrics - output as plot
  names(clustout_euclidean) <- paste0(2:10)
  cvi_out = t(sapply(clustout_euclidean, dtwclust::cvi, type = "internal"))
  cvi_out = data.frame(cvi_out) %>% mutate(n_clusters = as.numeric(names(clustout_euclidean) ))
  cvi_out = cvi_out %>% pivot_longer(!n_clusters, names_to = 'type', values_to = 'score')
  cvi_out %>% 
    ggplot() +
    geom_line(aes(x=n_clusters, y=score, color=type))+
    facet_wrap(~type, scale='free')
  
  ## Check max and mean cluster sizes at different n 
  plot(sapply(1:9, function(x){mean(clustout_euclidean[[x]]@clusinfo$size * clustout_euclidean[[x]]@clusinfo$av_dist/sum(clustout_euclidean[[x]]@clusinfo$size))})) 
  points(sapply(1:9, function(x){max(clustout_euclidean[[x]]@clusinfo$size * clustout_euclidean[[x]]@clusinfo$av_dist/sum(clustout_euclidean[[x]]@clusinfo$size))})) 
  
  ## output the dtw distances betweeen deviations to estimate cluster strengths
  distmat = data.table(clustout_euclidean[[7]]@distmat)
  aras = rownames(clustout_euclidean[[7]]@distmat)
  distmat[, ara_1 := aras]
  distances = melt(distmat, id.vars = 'ara_1', variable.name = 'ara_2', value.name = 'Dist')
  
  clusters = data.table(ara = aras, cluster = clustout_euclidean[[7]]@cluster)
  
  distances_clust = merge(distances, clusters, by.x = 'ara_1', by.y = 'ara')
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
  cent_aras=names(clustout_euclidean[[7]]@centroids)
  # merged to spatial data and plot for records 
  tt_cluster_sf = merge(shapes$ara_shp, unique(tt_cluster[,c('ara_id', 'ara_name', 'cluster', 'strong_clust')], by=c('ara_id')))
  
  
  map= 
    ggplot() + 
    geom_sf(data=tt_cluster_sf, aes(fill = as.character(cluster)))+ 
    geom_sf(data=tt_cluster_sf %>% filter(ara_name %in% cent_aras), color='limegreen', linewidth=1, fill=NA)# + facet_wrap(~period)
  
  
  
  ggsave(paste0('../plots/map_period_', peri, '.png'))
  
 
  map= 
    ggplot() + 
    geom_sf(data=tt_cluster_sf %>% filter(strong_clust==T), aes(fill = as.character(cluster))) + 
    geom_sf(data=tt_cluster_sf %>% filter(strong_clust==F), aes(fill = as.character(cluster)), alpha=0.2) + 
    geom_sf(data=tt_cluster_sf %>% filter(ara_name %in% cent_aras), color='limegreen', linewidth=1, fill=NA)# + facet_wrap(~period)
  
  
  ggsave(paste0('../plots/map_period_strong_', peri, '.png'))
  
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
  
  
  
  ggsave(paste0('../plots/series_period_', peri, '.png'))
  
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
  
  
  
  ggsave(paste0('../plots/series_period_strong_', peri, '.png'))
  
  
  av_dists = data.table(clustout_euclidean[[7]]@clusinfo) 
  colnames(av_dists) = c('No. in cluster', 'Average distance')
  av_dists[, cluster := 1:nrow(av_dists)]
  av_dists[, period := peri]
  av_dists_all = rbind(av_dists_all, av_dists)
  
  tt_cluster[, period:=peri]
  tt_cluster_all = rbind(tt_cluster_all, tt_cluster)
  centroid_distances = data.table()
  for(n in cent_aras){
    for(m in cent_aras){
      centroid_distances_pair = data.table(
        ara_1 = c(n), 
        ara_2 = c(m), 
        cluster_1 = clusters[ara==n]$cluster,
        cluster_2 = clusters[ara==m]$cluster,
        distance = c(dtwclust::dtw2(tt_list[n], tt_list[m]))$distance, 
        period = c(peri))
      centroid_distances = rbind(centroid_distances, centroid_distances_pair)
      
    }
    centroid_distances_all = rbind(centroid_distances_all, centroid_distances)
  }
}

tt_cluster_all
av_dists_all %>% ggplot() + 
  ggbeeswarm::geom_beeswarm(aes(x =period, y=`Average distance`, color=period))

centroid_distances_all[!(ara_1 == ara_2)] %>% ggplot() + 
  geom_point(aes(x = rank(distance), y=`distance`, color=as.character(period)))


centroid_distances_all_av = merge(centroid_distances_all, av_dists_all[,c('Average distance', 'cluster', 'period')], by.x = c('cluster_1', 'period'), by.y = c('cluster', 'period'))


centroid_distances_all_av[, rat_av_cen_dist := distance/`Average distance`]


centroid_distances_all_av[!(ara_1 == ara_2)] %>% ggplot() + 
  geom_point(aes(x = rank(rat_av_cen_dist), y=rat_av_cen_dist, color=as.character(period)))


centroid_distances_all_av[!(ara_1 == ara_2) & !is.infinite(rat_av_cen_dist)] %>% ggplot() + 
  geom_violin(aes(x = period, y=rat_av_cen_dist, color=as.character(period))) + 
  geom_point(aes(x = period, y=rat_av_cen_dist) )


write.csv(centroid_distances_all_av, '../outputs/dtw_results/centroid_distances_all_av.csv')

write.csv(av_dists_all, '../outputs/dtw_results/av_dists_all.csv')

write.csv(tt_cluster_all, '../outputs/dtw_results/tt_cluster_all.csv')




tt_cluster_sf_all = merge(shapes$ara_shp, unique(tt_cluster_all[,c('ara_id', 'ara_name', 'cluster', 'strong_clust', 'period')]), by=c('ara_id'))
map_all= 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='darkgray', linewidth=0.5, alpha=0.5) +
  geom_sf(data=tt_cluster_sf_all, aes(fill = as.character(cluster)), color=NA, alpha=0.8) + 
  facet_wrap(~period, nrow=1) + 
  geom_sf(data=shapes$see_shp, fill='midnightblue', color=NA)+ 
  theme_minimal() + 
  scale_fill_discrete(name='Cluster') + 
  ggtitle('A')

ggsave(paste0('../plots/map_all.png'))

map_all_strong= 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='darkgray', linewidth=0.5, alpha=0.5) +
  geom_sf(data=tt_cluster_sf_all %>% filter(strong_clust==T), aes(fill = as.character(cluster)), color=NA, alpha=0.8) + 
  geom_sf(data=tt_cluster_sf_all %>% filter(strong_clust==F), aes(fill = as.character(cluster)), color=NA, alpha=0.2) + 
  facet_wrap(~period, nrow=1) + 
  geom_sf(data=shapes$see_shp, fill='midnightblue', color=NA)+ 
  theme_minimal() + 
  scale_fill_discrete(name='Cluster') + 
  ggtitle('B')


ggsave(paste0('../plots/map_all_strong.png'))


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


unique_distances = unique(centroid_distances_all[!(cluster_1 == cluster_2), c('cluster_1', 'cluster_2', 'distance', 'period')]) 


unique_distances[cluster_2 > cluster_1, distance:=NA]

cluster_distance_mat = unique_distances[!is.na(distance)] %>% ggplot() + 
  geom_tile(aes(x=cluster_1, y=cluster_2, fill=distance), )+
  geom_text(aes(x=cluster_1, y=cluster_2,label = round(distance,3)), color = "white", size = 4) +
  scale_fill_viridis_c(trans='reverse')+
  facet_wrap(~period, nrow=1) + 
  ggtitle('C') + 
  theme_minimal()


joint_plot = map_all/cluster_distance_mat + plot_layout(heights=c(5,2))

ggsave(paste0('../plots/map_all_distances.png'), width=20, height=10, plot = joint_plot)

joint_plot_sc = map_all_strong/cluster_distance_mat + plot_layout(heights=c(5,2))

ggsave(paste0('../plots/map_all_distances_sc.png'), width=20, height=10, plot = joint_plot_sc)

joint_plot_all = map_all/map_all_strong/cluster_distance_mat + plot_layout(heights=c(1,1,1))
ggsave(paste0('../plots/map_all_distances_all.png'), width=20, height=10, plot = joint_plot_all)



g_clust = tt_cluster_all %>% filter(period==1 & strong_clust) %>% 
  ggplot() +
  geom_hline(yintercept=1,linetype=2,alpha=.5) +
  #geom_ribbon(aes(x=date,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
  geom_line(aes(x=date,y=exp(norm_mean), color=NUTS2_name, group=ara_name), alpha=0.5) +
  #geom_line( data = tt_cluster_all %>% filter(period==1 & ara_name %in% cent_aras), aes(x=date,y=exp(norm_mean), group=ara_name), color='limegreen', linetype=2, linewidth=2) +
  facet_wrap(~cluster) +
  scale_colour_discrete() +
  #scale_fill_discrete(guide="none") +
  scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
  coord_cartesian(ylim=c(.05,20)) +
  labs(title=paste0("Clusters"),x="Day",y="Relative viral load by ARA") + 
  theme(axis.text.x = element_text(angle=45,hjust=1))
