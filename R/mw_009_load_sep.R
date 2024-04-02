#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load SEP data
# creation: jriou
# init date: 2023-10-20
#:::::::::::::::::::::::::::::

mw_009_load_sep = function(shp) {
  # compute SEP weights
  swissboundaries_BFS_NATION = st_read('data/spatial/bfs/swissboundaries3d_2023-01_2056_5728.shp/swissBOUNDARIES3D_1_4_TLM_LANDESGEBIET.shp')
  bfs_crs = st_crs(swissboundaries_BFS_NATION)
  pop_data = fread('data/population_statistics/STATPOP2021.csv')
  pop_data_slim = pop_data[,c('E_KOORD', 'N_KOORD',  'B21BTOT')]
  
  pop_sf = st_as_sf(pop_data_slim, coords = c("E_KOORD", "N_KOORD"), 
                    crs = bfs_crs, remove=FALSE) 
  
  sep_sf = readRDS("../../02_data/sep3/ssep3_user_geo.Rds") %>% 
    sf::st_transform(crs=bfs_crs)
  
  joint_sf <- st_join(sep_sf, pop_sf, join = st_nearest_feature)
  
  joint_dt = data.table(st_drop_geometry(joint_sf))
  # head(joint_dt)
  joint_dt[, hect_id := paste0(E_KOORD, N_KOORD)]
  joint_dt[, sep_pop_weight := mean(B21BTOT)/.N, by=c('hect_id')]
  
  sep_weights = joint_dt[, c('gisid', 'sep_pop_weight')]
  
  # load and transform to WGS 84
  sep = readRDS("../../02_data/sep3/ssep3_user_geo.Rds") %>% 
    sf::st_transform(crs=4326)
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
      dplyr::mutate(ara_n=as.character(i)) %>% 
      dplyr::relocate(ara_n)
    out = dplyr::bind_rows(out,tmp)
  }
  
  
  return(out)
}
  