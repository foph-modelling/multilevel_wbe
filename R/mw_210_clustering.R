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
  alldays = unique(tt$ID)
  ndays = length(alldays)
  
  tt = tt %>%  
    dplyr::bind_cols(day=rep(alldays,nara),
                     ara1=rep(1:nara,each=ndays)) %>% 
    dplyr::left_join(corr,by = join_by(ara1)) %>% 
    dplyr::left_join(corr_days,by="day") %>% 
    as_tibble()

  tt2 = dplyr::filter(tt,ara_name %in% extremes$ara_name) %>% 
    dplyr::mutate(ara_name=factor(ara_name,levels=extremes$ara_name)) 

g1 = tt %>% filter(date<lubridate::ymd(20230101) )%>% 
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


write.csv(tt, '../data/true_shapes.csv')

tt_wide = data.table(dcast(tt, formula = day ~ ara_name, value.var = 'mean'))

tt_list = as.list(tt_wide[,-c('day')])



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

cluster_df[, cluster_tad:=clustout_tadpole@cluster]

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


