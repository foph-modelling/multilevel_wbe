#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for relative VL by ARA
# creation: jriou 
# init date: 2023-06-14
#:::::::::::::::::::::::::::::

mw_130_map_relative_vl = function(model,corr,shp) {
  
  # model = ma5.3.1
  # corr = corr_all_ara
  # shp = shapes
  
  # extract ara-level intercepts
  tt = model$summary.random$ara1 %>% 
    as_tibble() %>% 
    dplyr::transmute(ara1=ID,
                     `exp(beta)`=round(exp(mean),2),
                     `0.025quant`=round(exp(`0.025quant`),2),
                     `0.975quant`=as.numeric(format(round(exp(`0.975quant`),2),scientific=FALSE))) %>% 
    dplyr::left_join(corr,by = join_by(ara1)) %>% 
    dplyr::arrange(`exp(beta)`) %>% 
    dplyr::mutate(rank=as.factor(row_number()))
  
  # merge with map
  mm = shp$ara_shp %>% 
    dplyr::left_join(tt, by = join_by(ara_id)) %>% 
    dplyr::filter(!is.na(`exp(beta)`))
  
  # plot deviations
  g1 = ggplot(tt) +
    geom_pointrange(aes(x=rank,y=`exp(beta)`,ymin=`0.025quant`,ymax=`0.975quant`,colour=NUTS2_name)) +
    scale_x_discrete(labels=tt$ara_name) +
    scale_y_continuous(trans="pseudo_log",breaks=c(0,1,2,3,4,5,10)) +
    geom_hline(yintercept=1,linetype=2) +
    coord_flip() +
    theme(legend.position=c(.8,.2),
          legend.background = element_blank(),
          axis.text.y=element_text(size=6)) +
    labs(x=NULL,colour="NUTS2 region",title="Relative viral load by ARA compared to average") 
  
  # plot map
  g2 = ggplot() +
    geom_sf(data=shp$canton_shp,fill="grey95",colour="grey70") +
    geom_sf(data=shp$see_shp,fill="white") +
    geom_sf(data=mm,colour="black",aes(fill=`exp(beta)`)) +
    scale_fill_viridis_c(trans="log10",limits=c(0.1,10)) +
    labs(title="Relative viral load by ARA compared to average",fill="Relative VL") 
  
  return(g2)
}