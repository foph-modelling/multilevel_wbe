#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for missing and available wastewater data
# creation: jriou 
# init date: 2023-05-09
#:::::::::::::::::::::::::::::

mw_110_map_missing = function(ww,shp,identify_eawag=FALSE) {
  # gather info
  tt = ww %>% 
    dplyr::group_by(ara_id,NUTS2_name) %>% 
    dplyr::count() 
  
  # merge
  tt = shp$ara_shp %>% 
    dplyr::left_join(tt,by = join_by(ara_id)) %>% 
    filter(!is.na(n)) %>% 
    mutate(n_bin=case_when(n<200 ~ "<200",
                           n>=200 & n<500 ~ "200-500",
                           n>=500 ~ "500+"))
  
  # plot map
  # g = ggplot() +
  #   geom_sf(data=shp$canton_shp,fill="grey95",colour="grey70") +
  #   geom_sf(data=shp$see_shp,fill="white") +
  #   geom_sf(data=tt,colour="black",aes(fill=n)) +
  #   scale_fill_viridis_c(breaks=c(100,300,500,700),option="magma") +
  #   labs(fill="Number of measurements")   +
  #   theme(legend.position="bottom")
  # 
  # labels = tibble(NUTS2_name=c("Lake Geneva"),
  #               x=6,
  #               y=46)
  
  g = ggplot() +
    geom_sf(data=shp$canton_shp,fill="white",colour="black") +
    geom_sf(data=tt,aes(fill=n_bin),colour=NA) +
    geom_sf(data=tt,fill=NA,aes(colour=NUTS2_name),linewidth=.8) +
    geom_sf(data=filter(tt,NUTS2_name=="Central"),fill=NA,aes(colour=NUTS2_name),linewidth=.8) +
    # scale_fill_distiller(palette="Greys",direction=1,breaks=c(200,400,600)) +
    scale_fill_manual(values=c("grey80","grey50","grey10"),guide="none") +
    labs(colour=NULL,fill="Measurements")   +
    theme(legend.position="bottom")

  return(g)
}