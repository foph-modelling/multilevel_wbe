#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load shape files
# creation: jriou & mwagner
# init date: 2023-05-09
#:::::::::::::::::::::::::::::

mw_005_load_shp = function() {
  # load ----
  ara_shp = sf::read_sf("data/spatial/230214_ARA_BAG/230214_ARA_BAG.shp") %>% 
    dplyr::select(ara_id=ara_id,geometry)
  canton_shp = sf::read_sf("data/spatial/shp.shp")
  see_shp = readRDS("data/spatial/se.Rds")
  
  res = list(ara_shp=ara_shp,canton_shp=canton_shp,see_shp=see_shp)
  
  return(res)
  
  if(FALSE) {
    ggplot() +
      geom_sf(data=res$canton_shp,fill="grey95") +
      geom_sf(data=res$see_shp,fill="skyblue") +
      geom_sf(data=res$ara_shp,colour="red",fill="grey70")
  }
}