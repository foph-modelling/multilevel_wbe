#' ---
#' title: "WBE for SARS-CoV-2 in Switzerland: Clustering viral load time series"
#' author: "James Munday"
#' date: "`r Sys.Date()`"
#' params:
#'    controls: controls
#' output:
#'    html_document:
#'      code_folding : hide
#'      number_sections: false
#'      highlight: pygments
#'      theme: cosmo
#'      link-citations: true
#' ---
#' ## Dendograms of Dynamic Time-Warping distance
#' 
#'
#' 
 for (peri in 1:5){
   clust = readRDS(paste0('../outputs/cluster_', peri, '.rds'))
   dend = as.dendrogram(clust[[9]])
   
   nuts = unique(data.table(tt1)[ara_name %in% clust[[9]]$labels, c('ara_name', 'NUTS2_name')][order(ara_name),])$NUTS2_name
   colors_to_use <- as.numeric(as.factor(nuts))
   colors_to_use
   # But sort them based on their order in dend:
   colors_to_use <- colors_to_use[order.dendrogram(dend)]
   colors_to_use
   
   #pdf(paste0('../plots/dendogram_',peri,'.pdf'), width = 20, height=10 )
   # Now we can use them
   dendextend::labels_colors(dend) <- colors_to_use
   #dendextend::leaves_colors(dend) <- colors_to_use
   plot(dend)
   
print(paste0("*Figure ", peri, ".* Dendogram of period ", peri))
   

   }
#' 
#' 
#' ## Maps of clusters for each period and each k
#' 
#' 
#' 
#' 

 #tt_cluster_sf_all = readRDS(file = '../outputs/tt_cluster_sf_all.rds')
 #
 #map_all= 
 #  ggplot() + 
 #  geom_sf(data=shapes$canton_shp, fill=NA, color='darkgray', linewidth=0.5, alpha=0.5) +
 #  geom_sf(data=tt_cluster_sf_all, aes(fill = as.character(cluster)), color=NA, alpha=0.8) + 
 #  facet_grid(nclust~period) + 
 #  geom_sf(data=shapes$see_shp, fill='midnightblue', color=NA)+ 
 #  theme_minimal() + 
 #  scale_fill_discrete(name='Cluster') + 
 #  ggtitle('A') +
 #  theme(
 #    panel.grid = element_blank(),
 #    axis.text = element_blank(), 
 #    axis.ticks = element_blank())
 #
 #
 #map_all
 #
 #
 #g_clust = tt_cluster_all %>% filter((period==1 & nclust==8) |  (period==2 & nclust==5) | (period %in% c(3,4,5) & nclust==3)) %>% 
 #  ggplot() +
 #  geom_hline(yintercept=1,linetype=2,alpha=.5) +
 #  #geom_ribbon(aes(x=date,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
 #  geom_line(aes(x=date,y=exp(norm_mean), color=NUTS2_name, group=ara_name), alpha=0.5) +
 #  #geom_line( data = tt_cluster_all %>% filter(period==1 & ara_name %in% cent_aras), aes(x=date,y=exp(norm_mean), group=ara_name), color='limegreen', linetype=2, linewidth=2) +
 #  facet_grid(cluster~period, scales = 'free') +
 #  scale_colour_discrete() +
 #  #scale_fill_discrete(guide="none") +
 #  scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
 #  coord_cartesian(ylim=c(.05,20)) +
 #  labs(x="Day",y="Relative viral load by ARA") + 
 #  theme(axis.text.x = element_text(angle=45,hjust=1)) +
 #  theme_minimal()
 #
 #g_clust
 #
 #