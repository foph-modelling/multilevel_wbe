#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for missing and available wastewater data
# creation: jriou 
# init date: 2023-05-09
#:::::::::::::::::::::::::::::

mw_110_map_missing = function(ww,shp,identify_eawag=FALSE) {
  # gather info
  tt = ww %>% 
    dplyr::filter(measurement==1) %>% 
    dplyr::group_by(ara_id) %>% 
    dplyr::count() 
  
  # merge
  tt = shp$ara_shp %>% 
    dplyr::left_join(tt,by = join_by(ara_id)) %>% 
    filter(!is.na(n))
  
  # plot map
  g = ggplot() +
    geom_sf(data=shp$canton_shp,fill="grey95",colour="grey70") +
    geom_sf(data=shp$see_shp,fill="skyblue") +
    geom_sf(data=tt,colour="black",aes(fill=n)) +
    scale_fill_viridis_c() +
    labs(fill="Measurements")

  return(g)
}