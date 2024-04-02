#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load SEP data
# creation: jriou
# init date: 2023-10-20
#:::::::::::::::::::::::::::::

mw_009_load_sep = function(shp) {
  # load and transform to WGS 84
  sep = readRDS("../../02_data/sep3/ssep3_user_geo.Rds") %>% 
    sf::st_transform(crs=4326)
  sep_weights = fread('data/sep_weights.csv')
  sep = sep %>% mutate(pop_weight = sep_weights$sep_pop_weight)
  
  # find SEP by ARA
  jj = shp$ara_shp %>% 
    sf::st_contains(sep) 
  
  # global sd
  glob_sd = sqrt(Hmisc::wtd.var(sep$ssep3, sep$pop_weight))

  # find average SEP by ARA
  out = NULL
  n_ara = nrow(shp$ara_shp)
  sep2 = as_tibble(sep)
  for(i in 1:n_ara) {
    tmp = dplyr::filter(sep2,gisid %in% jj[[i]]) %>% 
      dplyr::summarise(ssep3_med=Hmisc::wtd.quantile(ssep3, pop_weight, probs = 0.5)[[1]],
                       ssep3_mean=Hmisc::wtd.mean(ssep3, pop_weight),
                       ssep3_sd=sqrt(Hmisc::wtd.var(ssep3, pop_weight)),
                       ssep3_icc=ssep3_sd/glob_sd,
                       ssep3_min=min(ssep3),
                       ssep3_max=max(ssep3)) %>% 
      dplyr::mutate(ara_id=as.character(i)) %>% 
      dplyr::relocate(ara_id)
    out = dplyr::bind_rows(out,tmp)
  }
  
  
  return(out)
}
  