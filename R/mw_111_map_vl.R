#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for SARS-CoV-2 viral load in wastewater data
# creation: jriou 
# init date: 2023-05-09
#:::::::::::::::::::::::::::::

mw_111_map_vl = function(ww,shp) {
  # gather info
  tt = ww %>% 
    dplyr::filter(below_lod==0,below_loq==0) %>% 
    dplyr::group_by(ara_id,period,week) %>% 
    dplyr::summarise(vl=mean(vl,na.rm=TRUE),.groups="drop") %>% 
    dplyr::group_by(ara_id,period) %>% 
    dplyr::summarise(vl=mean(vl),.groups="drop")
  
  # merge
  tt = shp$ara_shp %>% 
    dplyr::left_join(tt, by = join_by(ara_id)) %>% 
    filter( !is.na(period))
  
  # range
  vl_range = ww %>% 
    dplyr::group_by(ara_id,period,week) %>% 
    dplyr::summarise(vl=mean(vl,na.rm=TRUE),.groups="drop") %>% 
    dplyr::filter(!is.na(vl)) %>% 
    dplyr::summarise(min=min(vl,na.rm=TRUE),max=max(vl,na.rm=TRUE))

  # plot map
  g = ggplot() +
    geom_sf(data=shp$canton_shp,fill="grey95",colour="grey70") +
    geom_sf(data=shp$see_shp,fill="white") +
    geom_sf(data=tt,colour="black",aes(fill=vl)) +
    scale_fill_viridis_c(trans="log10") +
    facet_wrap(~period)
  
  return(g)
}