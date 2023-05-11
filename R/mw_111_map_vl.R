#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for SARS-CoV-2 viral load in wastewater data
# creation: jriou 
# init date: 2023-05-09
#:::::::::::::::::::::::::::::

mw_111_map_vl = function(ww,shp,lower_limit=0) {
  # gather info
  tt = ww %>% 
    dplyr::group_by(ara_id,period,week) %>% 
    dplyr::summarise(vl=mean(vl,na.rm=TRUE),.groups="drop") %>% view()
    dplyr::filter(!is.na(vl),vl>lower_limit) %>% 
    dplyr::group_by(ara_id,period) %>% 
    dplyr::summarise(vl=max(vl),.groups="drop")
  
  # merge
  tt = shp$ara_shp %>% 
    dplyr::left_join(tt, by = join_by(ara_id)) %>% 
    filter(!is.na(n), !is.na(period))
  
  # plot map
  g = ggplot() +
    geom_sf(data=shp$canton_shp,fill="grey95",colour="grey70") +
    geom_sf(data=shp$see_shp,fill="skyblue") +
    geom_sf(data=tt,colour="black",aes(fill=vl)) +
    scale_fill_viridis_c() +
    facet_wrap(~period)
  
  return(g)
}