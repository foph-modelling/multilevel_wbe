library(patchwork)
library(tidyverse)
library(data.table)

model = readRDS(paste0("../",controls$savepoint,"ma5.3.2.rds"))

ww1 = readRDS(file=paste0("../",controls$savepoint,"ww_all.rds"))

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
tt1 = model$summary.random$day1
alldays = unique(tt1$ID)
ndays = length(alldays)

tt1 = tt1 %>%  
  dplyr::bind_cols(day=rep(alldays,nara),
                   ara1=rep(1:nara,each=ndays)) %>% 
  dplyr::left_join(corr,by = join_by(ara1)) %>% 
  dplyr::left_join(corr_days,by="day") %>% 
  as_tibble()

  #tt2 = dplyr::filter(tt1,ara_name %in% extremes$ara_name) %>% 
  #  dplyr::mutate(ara_name=factor(ara_name,levels=extremes$ara_name)) 
  
  
tt_clusters = tibble()
  for(peri in 1:4){
    tt = tt1  %>% filter(period %in% c(3))
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
                                           distance = 'dtw2',
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
    
    cluster_df[, cluster_tad:=clustout_tadpole[[4]]@cluster]
    
    tt_cluster = merge(tt, cluster_df, by=c('ara_name'))
    
    dtwclust::interactive_clustering(series = tt1)
    
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
    
    tt_cluster %>% mutate(period=peri)
    tt_clusters = rbind(tt_clusters, tt_cluster)
  }

tt_cluster_sf_tad = merge(shapes$ara_shp, unique(tt_clusters[,c('ara_id', 'ara_name', 'cluster_tad', 'period')], by=c('ara_id')))
tt_cluster_sf_tad = merge(shapes$ara_shp, unique(tt_cluster[,c('ara_id', 'ara_name', 'cluster_tad')], by=c('ara_id')))

map= tt_cluster_sf_tad %>% ggplot() + geom_sf(aes(fill = as.character(cluster_tad)))# + facet_wrap(~period)




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


clustout_euclidean = dtwclust::tsclust(series=dtwclust::zscore(tt_list), 
                                       type='partitional', 
                                       distance = 'dtw2',
                                       k=10, 
                                       return.objects=TRUE)



cluster_df[, cluster_euc:=clustout_euclidean@cluster]


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


clustout_tadpole = dtwclust::tsclust(series=dtwclust::zscore(tt_list), 
                                     type='tadpole',
                                     k=2:25, 
                                     control = dtwclust::tadpole_control(dc=6, window.size = 18L),
                                     return.objects=TRUE)

names(clustout_tadpole) <- paste0(2:25)
cvi_out = t(sapply(clustout_tadpole, dtwclust::cvi, type = "internal"))
cvi_out = data.frame(cvi_out) %>% mutate(n_clusters = as.numeric(names(clustout_tadpole) ))
cvi_out = cvi_out %>% pivot_longer(!n_clusters, names_to = 'type', values_to = 'score')
cvi_out %>% 
  ggplot() +
  geom_line(aes(x=n_clusters, y=score, color=type))+
  facet_wrap(~type, scale='free')



clustout_tadpole = dtwclust::tsclust(series=dtwclust::zscore(tt_list), 
                                     type='tadpole',
                                     k=10, 
                                     control = dtwclust::tadpole_control(dc=8, window.size =10L),
                                     return.objects=TRUE)

cluster_df = data.table(
  ara_name = names(tt_list)
)

cluster_df[, cluster_tad:=clustout_euclidean@cluster]

tt_cluster = merge(tt, cluster_df, by=c('ara_name'))


g_clust = tt_cluster %>% filter(date<lubridate::ymd(20230101) )%>% 
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

g_clust



tt_cluster_sf_tad = merge(shapes$ara_shp, unique(tt_cluster[,c('ara_id', 'ara_name', 'cluster_tad')], by=c('ara_id')))

map= tt_cluster_sf_tad %>% ggplot() + geom_sf(aes(fill = as.character(cluster_tad)))

g_clust + map

tt_cluster_sf_spat = merge(shapes$ara_shp, unique(tt_cluster[,c('ara_id', 'cluster_spat')], by=c('ara_id')))

map= tt_cluster_sf_spat %>% ggplot() + geom_sf(aes(fill = as.character(cluster_spat)))

tt_cluster_sf_te = merge(shapes$ara_shp, unique(tt_cluster[,c('ara_id', 'cluste')], by=c('ara_id')))

map= tt_cluster_sf_te %>% ggplot() + geom_sf(aes(fill = as.character(cluste)))

tt_cluster_sf_euc = merge(shapes$ara_shp, unique(tt_cluster[,c('ara_id', 'cluster_euc')], by=c('ara_id')))
#
map = tt_cluster_sf_euc %>% ggplot() + geom_sf(aes(fill = as.character(cluster_euc)))
#

plot_list = list()
for(clust in 1:10){ 
  
  p = 
    ggplot() + 
    geom_sf(data = tt_cluster_sf_spat, fill='lightgray') + 
    geom_sf(data = tt_cluster_sf_spat %>% filter(cluster_spat == clust), aes(fill = as.character(cluster_spat))) + 
    theme(legend.position = 'none')
  plot_list[[clust]] = p
  
  }
patchwork::wrap_plots(plot_list, ncol=1)

### Transfer Entropy

tt_ind = 1:length(tt_list)

te = RTransferEntropy::transfer_entropy(tt_list[[6]], tt_list[[7]], lx=1, ly=1)

future::plan(future::multisession, workers = 6)
te_mat_list = purrr::map2(expand.grid(tt_ind, tt_ind)[,1], expand.grid(tt_ind, tt_ind)[,2], function(x,y){sum(RTransferEntropy::transfer_entropy(tt_list[[x]], tt_list[[y]])$coef[,1])})

te_dt = rbind(tt_ind, unlist())

image(matrix(unlist(te_mat_list), nrow=length(tt_list)))









## ------  compare distances in clusters ----------

av_dists_all = data.table()
tt_cluster_all = data.table()
centroid_distances_all = data.table()

for (peri in 1:4){
  tt = data.table(tt1  %>% filter(period %in% c(peri)))
  tt_masked = left_join(tt, ww %>% select(date, ara_name, vl), by=c('date', 'ara_name'))
  mask_in_period = tt_masked %>% 
    group_by(ara_id) %>%
    summarise(x = sum(vl, na.rm = TRUE))
  
  mask_dates= tt_masked %>% 
    group_by(ara_id, day) %>%
    mutate(include = !is.na(vl))
  
  
  mask_from_na = mask_dates %>% 
    group_by(ara_id) %>%
    mutate(x = sum(include, na.rm = TRUE)/length(include))
  
  aras_to_include_na = unique(mask_from_na %>% filter(x > 0.3) %>% select(ara_id))
  
  
  aras_to_include = mask_in_period %>% filter(x != 0) %>% select(ara_id)
  
  tt =  tt %>% filter((ara_id %in% unlist(aras_to_include)))
  tt =  tt %>% filter((ara_id %in% unlist(aras_to_include_na)))
  
  
  tt_wide = data.table(dcast(tt, formula = day ~ ara_name, value.var = 'mean'))
  
  tt_list = as.list(tt_wide[,-c('day')])
  
  
  
  
  clustout_euclidean = dtwclust::tsclust(series=dtwclust::zscore(tt_list), 
                                         type='hierarchical', 
                                         distance = 'dtw2',
                                         k=2:25, 
                                         return.objects=TRUE)
  names(clustout_euclidean) <- paste0(2:25)
  cvi_out = t(sapply(clustout_euclidean, dtwclust::cvi, type = "internal"))
  cvi_out = data.frame(cvi_out) %>% mutate(n_clusters = as.numeric(names(clustout_euclidean) ))
  cvi_out = cvi_out %>% pivot_longer(!n_clusters, names_to = 'type', values_to = 'score')
  cvi_out %>% 
    ggplot() +
    geom_line(aes(x=n_clusters, y=score, color=type))+
    facet_wrap(~type, scale='free')
  
  
  plot(sapply(1:24, function(x){mean(clustout_euclidean[[x]]@clusinfo$size * clustout_euclidean[[x]]@clusinfo$av_dist/sum(clustout_euclidean[[x]]@clusinfo$size))})) 
  points(sapply(1:24, function(x){max(clustout_euclidean[[x]]@clusinfo$size * clustout_euclidean[[x]]@clusinfo$av_dist/sum(clustout_euclidean[[x]]@clusinfo$size))})) 
  
  
  distmat = data.table(clustout_euclidean[[7]]@distmat)
  aras = rownames(clustout_euclidean[[7]]@distmat)
  distmat[, ara_1 := aras]
  distances = melt(distmat, id.vars = 'ara_1', variable.name = 'ara_2', value.name = 'Dist')
  
  clusters = data.table(ara = aras, cluster = clustout_euclidean[[7]]@cluster)
  
  distances_clust = merge(distances, clusters, by.x = 'ara_1', by.y = 'ara')
  distances_clust = distances_clust[, cluster_1 := cluster][,-'cluster']
  distances_clust = merge(distances_clust, clusters, by.x = 'ara_2', by.y = 'ara')
  distances_clust = distances_clust[, cluster_2 := cluster][,-'cluster']
  
  
  compare_distances = function(distances_clust, clust, ara){
    distances_to_common_clust = distances_clust[ara_1 == ara & ara_2 != ara & cluster_2 == clust, ]$Dist
    distances_to_other_clusts = distances_clust[ara_1 == ara  & cluster_2 != clust, ]$Dist
    quantile(distances_to_common_clust, probs = 0.5) < quantile(distances_to_other_clusts, probs = 0.05)
  }
  
  
  
  clusters[, strong_clust := mapply(x=clusters$ara, y=clusters$cluster, FUN = function(x,y){compare_distances(distances_clust,y,x)}) ]
  
  
  tt_cluster = merge(tt, clusters, by.x=c('ara_name'), by.y='ara')
  
  
  
  
  
  tt_cluster_sf = merge(shapes$ara_shp, unique(tt_cluster[,c('ara_id', 'ara_name', 'cluster', 'strong_clust')], by=c('ara_id')))
  map= 
    ggplot() + 
    geom_sf(data=tt_cluster_sf, aes(fill = as.character(cluster)))+ 
    geom_sf(data=tt_cluster_sf %>% filter(ara_name %in% cent_aras), color='limegreen', linewidth=1, fill=NA)# + facet_wrap(~period)
  
  
  
  ggsave(paste0('../plots/map_period_', peri, '.png'))
  
  cent_aras=names(clustout_euclidean[[7]]@centroids)
  map= 
    ggplot() + 
    geom_sf(data=tt_cluster_sf %>% filter(strong_clust==T), aes(fill = as.character(cluster))) + 
    geom_sf(data=tt_cluster_sf %>% filter(strong_clust==F), aes(fill = as.character(cluster)), alpha=0.2) + 
    geom_sf(data=tt_cluster_sf %>% filter(ara_name %in% cent_aras), color='limegreen', linewidth=1, fill=NA)# + facet_wrap(~period)
  
  
  ggsave(paste0('../plots/map_period_strong_', peri, '.png'))
  
  tt_cluster[, norm_mean := dtwclust::zscore(mean), by=c('ara_id')]
  
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
  scale_fill_discrete(name='Cluster')

ggsave(paste0('../plots/map_all.png'))

map_all_strong= 
  ggplot() + 
  geom_sf(data=shapes$canton_shp, fill=NA, color='darkgray', linewidth=0.5, alpha=0.5) +
  geom_sf(data=tt_cluster_sf_all %>% filter(strong_clust==T), aes(fill = as.character(cluster)), color=NA, alpha=0.8) + 
  geom_sf(data=tt_cluster_sf_all %>% filter(strong_clust==F), aes(fill = as.character(cluster)), color=NA, alpha=0.2) + 
  facet_wrap(~period, nrow=1) + 
  geom_sf(data=shapes$see_shp, fill='midnightblue', color=NA)+ 
  theme_minimal() + 
  scale_fill_discrete(name='Cluster')


ggsave(paste0('../plots/map_all_strong.png'))


pal <- RColorBrewer::brewer.pal(8, name ="Set3")

for(peri in 1:4) {
  centroid_graph = igraph::graph_from_data_frame(unique(centroid_distances_all[period==peri & distance > 0,c('cluster_1', 'cluster_2', 'distance')]), directed = F)
  igraph::edge.attributes(centroid_graph)$weight <- log(1.0/(unique(centroid_distances_all[period==peri & distance > 0,c('cluster_1', 'cluster_2', 'distance')])$distance))
  igraph::E(centroid_graph)$weight
  plot(centroid_graph, 
       edge.arrow.size= 0, 
       edge.curved = TRUE, 
       edge.width = igraph::E(centroid_graph)$weight * 10,
       vertex.color=pal[ as.numeric(igraph::V(centroid_graph)$name)],
       main=paste0('period ' ,peri))
}

image(as.matrix(dcast(unique(centroid_distances_all[period==peri, c('cluster_1', 'cluster_2', 'distance')]), formula=cluster_1~cluster_2, value.var = 'distance')))


unique_distances = unique(centroid_distances_all[!(cluster_1 == cluster_2), c('cluster_1', 'cluster_2', 'distance', 'period')]) 


unique_distances[cluster_2 > cluster_1, distance:=NA]

cluster_distance_mat = unique_distances[!is.na(distance)] %>% ggplot() + 
  geom_tile(aes(x=cluster_1, y=cluster_2, fill=distance), )+
  geom_text(aes(x=cluster_1, y=cluster_2,label = round(distance,3)), color = "white", size = 4) +
  scale_fill_viridis_c(trans='reverse')+
  facet_wrap(~period, nrow=1) + 
  theme_minimal()


map_all/cluster_distance_mat + plot_layout(heights=c(5,2))

ggsave(paste0('../plots/map_all_distances.png'), width=20, height=10)






g_clust = tt_cluster_all %>% filter(period==1) %>% 
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
