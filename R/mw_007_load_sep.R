#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load SEP data
# creation: jriou
# init date: 2023-10-20
#:::::::::::::::::::::::::::::

mw_007_load_sep = function(shp) {
  # load and transform to WGS 84
  sep = readRDS("../../02_data/sep3/ssep3_user_geo.Rds") %>% 
    sf::st_transform(crs=4326)
  
  # find SEP by ARA
  jj = shapes$ara_shp %>% 
    sf::st_contains(sep) 

  # find average SEP by ARA
  n_ara = nrow(shapes$ara_shp)
  sep2 = as_tibble(sep)
  for(i in 1:n_ara) {
    dplyr::filter(sep2,gisid %in% jj[[i]])
  }
  
  do.call("rbind",jj) %>% 
    as_tibble()
  
  
}
  