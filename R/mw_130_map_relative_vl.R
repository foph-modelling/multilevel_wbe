#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for relative VL by ARA
# creation: jriou 
# init date: 2023-06-14
#:::::::::::::::::::::::::::::

mw_130_map_relative_vl = function(model,corr,shp,forestplot=FALSE, top=NULL,pprint=FALSE) {
  # model = ma5.3.1
  # corr = corr_all_ara
  # shp = shapes
  
  # nuts2 shapes
  
  shape_kt = shp$canton_shp %>%
    dplyr::mutate(kt=c("VD","VS","GE","BE","FR","SO","NE","JU","BS","BL","AG","ZH","GL",
                       "SH","AR","AI","SG","GR","TG","LU","UR","SZ","OW","NW","ZG","TI")) %>% 
    left_join(corr,by="kt") 
  shape_nuts = shape_kt %>% 
    group_by(NUTS2_name) %>% 
    summarise(geometry=sf::st_union(geometry))
  # shape_nuts = sf::read_sf("data/spatial/NUTS_RG_20M_2021_3035.shp/NUTS_RG_20M_2021_3035.shp") %>% 
    # filter(CNTR_CODE=="CH",LEVL_CODE==2)
  
  # extract ara-level intercepts
  tt = model$summary.random$ara1 %>% 
    as_tibble() %>% 
    dplyr::transmute(ara1=ID,
                     `exp(beta)`=round(exp(mean),2),
                     `0.025quant`=round(exp(`0.025quant`),2),
                     `0.975quant`=as.numeric(format(round(exp(`0.975quant`),2),scientific=FALSE))) %>% 
    filter(ara1<118) %>% 
    dplyr::left_join(corr,by = join_by(ara1)) %>% 
    dplyr::arrange(`exp(beta)`) %>% 
    dplyr::mutate(rank=as.factor(row_number()))
  
  if(!is.null(top)) {
    tt = dplyr::bind_rows(head(tt,top),tail(tt,top))
  }
  
  if(pprint) print(dplyr::bind_rows(head(tt,1),tail(tt,1)))
  # merge with map
  mm = shp$ara_shp %>% 
    dplyr::left_join(tt, by = join_by(ara_id)) %>% 
    dplyr::filter(!is.na(`exp(beta)`))
  
  # plot deviations
  g1 = ggplot(tt) +
    geom_pointrange(aes(x=rank,y=`exp(beta)`,ymin=`0.025quant`,ymax=`0.975quant`,colour=NUTS2_name)) +
    scale_x_discrete(labels=tt$ara_name) +
    scale_y_continuous(trans="pseudo_log",breaks=c(0,.5,1,2,3,4,5,10)) +
    scale_colour_discrete(guide="none") +
    geom_hline(yintercept=1,linetype=2) +
    coord_flip() +
    # theme(legend.position=c(.8,.2),
          # legend.background = element_blank(),
          # axis.text.y=element_text(size=6)) +
    labs(x=NULL,colour="NUTS2 region",title="Relative viral load by ARA compared to average") 
  
  # plot map
  g2 = ggplot() +
    geom_sf(data=shape_nuts,fill="grey95",colour="grey70") +
    geom_sf(data=shp$see_shp,fill="white") +
    geom_sf(data=mm,colour="black",aes(fill=`exp(beta)`)) +
    scale_fill_viridis_c(trans="log",breaks=c(.5,.7,1,1.5,2)) +
    labs(fill="Fold change in VL")  +
    theme(legend.position="bottom")
  
  mm %>% as.data.frame()
  # ggplot() +
  #   
  #   geom_sf(data=shape_nuts,fill=NA,colour="grey70") +
  #   geom_sf(data=shp$canton_shp,aes(fill=NAME),colour="grey70") +
  #   geom_sf(data=shp$see_shp,fill="white") +
  #   scale_fill_grey() +
  #   new_scale_fill() + 
  #   geom_sf(data=shp$see_shp,fill="white") +
  #   geom_sf(data=mm,colour="black",aes(fill=`exp(beta)`)) +
  #   # geom_sf(data=shape_nuts,colour="black",fill=NA,linewidth=1.2,alpha=.2) +
  #   scale_fill_viridis_c(trans="log",breaks=c(.5,.7,1,1.5,2)) +
  #   labs(fill="Fold change in VL")  +
  #   theme(legend.position="bottom")
  
  if(forestplot) {
    g = g1
  } else {
    g = g2
  }
  return(g)
}