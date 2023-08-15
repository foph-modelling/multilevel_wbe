#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load population from bfs hectare data and aggregate to plz and ARA
# creation: jmunday
# init date: 2023-08-15
#:::::::::::::::::::::::::::::

require(sf)
require(tidyverse)
require(data.table)

sf_use_s2(FALSE)

mw_007_load_hect_pop = function() {
  
  
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
  
  emp_data = data.table::fread('data/employment_statistics/STATENT_2020.csv')
  emp_data_slim = emp_data[,c('N_KOORD', 'E_KOORD', 'B08VZAT')]
  emp_data_slim_sf = st_as_sf(emp_data_slim, coords = c("E_KOORD", "N_KOORD"), 
                              crs = bfs_crs, remove=FALSE)
  
  st_transform(emp_data_slim_sf, crs = catchment_crs)
  
  plz_area = st_read('data/spatial/plz/PLZO_SHP_LV95/PLZO_PLZ.shp')
  plz_area = st_transform(plz_area, crs = catchment_crs)
  
  
  plz_catch_area = data.table()
  for (catch in catchments$ara_id){
    suppressMessages({
      plz_catch_area_part = data.table(st_filter(plz_area, catchments %>% filter(ara_id==catch), .predicate=st_intersects))
    })
    plz_catch_area_part[, ara_id:=catch]
    plz_catch_area = rbind(plz_catch_area, plz_catch_area_part)
  }
  
  plz_catch_area = st_as_sf(plz_catch_area, crs=st_crs(plz_area))
  
  
  
  # extract proportions of plzs actually in catchments - each plz/catchment combination is it's own polygon
  suppressWarnings({suppressMessages({
    intersections = st_intersection(catchments, plz_catch_area)
  })})
  
  
  # join population data and plz-catchment shapes data to extract all hectares within each polygon
  # NB the center of the hectare is used for location - this could be changed to include the proportion
  # of the hectare * population but this is quite involved.
  catch_pop_data = st_join(intersections, st_transform(pop_data_slim_sf, crs=catchment_crs))
  catch_pop_data_dt = data.table(catch_pop_data)[,c('ara_id', 'ara_name', 'PLZ', 'N_KOORD', 'E_KOORD', 'B21BTOT')]
  
  catch_pop_data_dt[, plz_catch_pop := sum(B21BTOT), by=c('ara_id', 'PLZ')]
  
  plz_catch_pops = unique(catch_pop_data_dt[, c('ara_id', 'ara_name', 'PLZ', 'plz_catch_pop')]) %>% filter(!is.na(plz_catch_pop))
  
  plz_catch_pops[, pop_total := sum(plz_catch_pop), by = c('ara_id')]
  
  ara_pop = unique(plz_catch_pops[, c('ara_id', 'pop_total')])
  
  return(ara_pop)
  
}
