
require(sf)
require(tidyverse)
require(data.table)

sf_use_s2(FALSE)

mw_012_calc_props_plz_catch = function() {
  
  
  # load shape objects
  ## load catchment geometry
  catchments = read_sf('data/spatial/230214_ARA_BAG/230214_ARA_BAG.shp')
  catchments = st_make_valid(catchments)
  catchment_crs = st_crs(catchments)
  
  ## load population data (Total population by hectare (2021) (100mx100m squares))
  swissboundaries_BFS_NATION = st_read('data/spatial/bfs/swissboundaries3d_2023-01_2056_5728.shp/swissBOUNDARIES3D_1_4_TLM_LANDESGEBIET.shp')
  bfs_crs = st_crs(swissboundaries_BFS_NATION)
  pop_data = fread('data/population_statistics/STATPOP2021.csv')
  pop_data_slim = pop_data[,c('E_KOORD', 'N_KOORD',  'B21BTOT')]
  
  pop_data_slim_sf = st_as_sf(pop_data_slim, coords = c("E_KOORD", "N_KOORD"), 
                              crs = bfs_crs, remove=FALSE) %>% st_transform(crs = catchment_crs)
  
  plz_area = st_read('data/spatial/plz/PLZO_SHP_LV95/PLZO_PLZ.shp')
  plz_area = st_transform(plz_area, crs = catchment_crs)
  
  
  # extract proportions of plzs actually in catchments - each plz/catchment combination is it's own polygon
  suppressWarnings({suppressMessages({
    intersections = st_intersection(catchments, plz_area)
  })})
  
  
  
  
  # join population data and plz-catchment shapes data to extract all hectares within each polygon
  # NB the center of the hectare is used for location - this could be changed to include the proportion
  # of the hectare * population but this is quite involved.
  catch_pop_data = st_join(st_transform(pop_data_slim_sf, crs=catchment_crs), intersections,)
  catch_pop_data_dt = data.table(catch_pop_data)[,c('ara_id', 'ara_name', 'PLZ', 'N_KOORD', 'E_KOORD', 'B21BTOT')]
  
  catch_pop_data_dt[, plz_catch_pop := sum(B21BTOT), by=c('ara_id', 'PLZ')]
  
  plz_catch_pops = unique(catch_pop_data_dt[, c('ara_id', 'ara_name', 'PLZ', 'plz_catch_pop')]) %>% filter(!is.na(plz_catch_pop))
  
  plz_catch_pops[, pop_total := sum(plz_catch_pop), by = c('ara_id')]
  
  suppressWarnings({suppressMessages({
    plz_pop_data = st_join(plz_area, st_transform(pop_data_slim_sf, crs=catchment_crs))
  })})
  
  plz_pop_data_dt = data.table(st_drop_geometry(plz_pop_data))
  
  plz_pop_data_dt[, plz_pop := sum(B21BTOT), by=c('PLZ')]
  
  plz_pops = unique(plz_pop_data_dt[, c('PLZ', 'plz_pop')])
  
  plz_catch_and_totals = merge(plz_catch_pops, plz_pops, how='left')
  plz_catch_and_totals[, prop_plz_in_catch:= plz_catch_pop/plz_pop]
  plz_catch_and_totals[, prop_catch_in_plz := plz_catch_pop/pop_total]
  
 
  
  return(plz_catch_and_totals)
  
}
