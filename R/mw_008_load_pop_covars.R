#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load population from bfs hectare data and aggregate to plz and ARA
# creation: jmunday
# init date: 2023-08-16
#:::::::::::::::::::::::::::::


mw_008_load_pop_covars = function(scale='ARA') {
  
  sf_use_s2(FALSE)
  
  # load shape objects
  ## load catchment geometry
  catchments = read_sf('data/spatial/230214_ARA_BAG/230214_ARA_BAG.shp')
  catchments = st_make_valid(catchments)
  catchment_crs = st_crs(catchments)
  
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
  
  
  ## load population data (Total population by hectare (2021) (100mx100m squares))
  swissboundaries_BFS_NATION = st_read('data/spatial/bfs/swissboundaries3d_2023-01_2056_5728.shp/swissBOUNDARIES3D_1_4_TLM_LANDESGEBIET.shp')
  bfs_crs = st_crs(swissboundaries_BFS_NATION)
  
  # extract proportions of plzs actually in catchments - each plz/catchment combination is it's own polygon
  suppressWarnings({suppressMessages({
    intersections = st_intersection(catchments, plz_area)
  })})
  
  ## load population data (Total population by hectare (2021) (100mx100m squares))
  swissboundaries_BFS_NATION = st_read('data/spatial/bfs/swissboundaries3d_2023-01_2056_5728.shp/swissBOUNDARIES3D_1_4_TLM_LANDESGEBIET.shp')
  bfs_crs = st_crs(swissboundaries_BFS_NATION)
  pop_data = fread('data/population_statistics/STATPOP2021.csv')
  
  # Define colulms of interest for model covariats 
  males = names(pop_data)[grepl('^B\\d{2}(BM)(0|1)', names(pop_data))] 
  fmale = names(pop_data)[grepl('^B\\d{2}(BW)(0|1)', names(pop_data))] 
  non_CH_EU = c('B21B28', 'B21B29', 'B21B30')
  select_cols = c('B21BTOT', males, fmale, non_CH_EU, "E_KOORD", "N_KOORD")
  
  
  # Trim data + make it an sf object to allow selection by polygon
  pop_data_slim = pop_data[,..select_cols]
  pop_data_slim_sf = st_as_sf(pop_data_slim, coords = c("E_KOORD", "N_KOORD"), 
                              crs = bfs_crs, remove=FALSE)
  
  
  ## Load employment data (Total FTE employment by hectare (2020) (100mx100m squares))
  emp_data = data.table::fread('data/employment_statistics/STATENT_2020.csv',fill = TRUE)
  
  # Define columns of interest for covariates
  
  
  emp_data_slim = emp_data[,c('N_KOORD', 'E_KOORD', 'B08VZAT')]
  emp_data_slim_sf = st_as_sf(emp_data_slim, coords = c("E_KOORD", "N_KOORD"), 
                              crs = bfs_crs, remove=FALSE)
  
  if(scale == 'ARA'){
    ### POPULATION DATA
    catch_pop_data = st_join(st_transform(pop_data_slim_sf, crs=st_crs(plz_catch_area)), intersections)
    catch_pop_data_no_geom = st_drop_geometry(catch_pop_data)[,c('ara_id', 'PLZ', select_cols )]
    catch_pop_data_dt = data.table(data.frame(catch_pop_data_no_geom))
    
    catch_cov_cols = c('ara_id', 'PLZ', 'B21BTOT', males, fmale, non_CH_EU)
    
    catch_pops_covs = catch_pop_data_dt[, ..catch_cov_cols]
    catch_pops_covs_summed = catch_pops_covs[, lapply(.SD, sum, na.rm=TRUE), by=c('ara_id'), .SDcols=c('B21BTOT', males, fmale, non_CH_EU) ] 
    
    
    
    under_20_cols = grep("^B\\d{2}(BM|BW)(0[1-3])$", c( males, fmale), value = TRUE)
    over_65_cols = grep("^B\\d{2}(BM|BW)(1[4-9]|2\\d|3[0-9])$",c(males, fmale), value = TRUE)
    
    
    catch_pops_covs_summed[, prop_under_20 := rowSums(.SD) / B21BTOT, .SDcols = under_20_cols  ]
    catch_pops_covs_summed[, prop_over_65 := rowSums(.SD) / B21BTOT, .SDcols = over_65_cols  ]
    
    catch_pops_covs_summed[, non_ch_eu := rowSums(.SD), .SDcols = non_CH_EU]
    catch_pops_covs_summed[, prop_non_ch_eu := non_ch_eu/B21BTOT]
    catch_pops_covs_summed$prop_non_ch_eu = nafill(catch_pops_covs_summed$prop_non_ch_eu,fill=0)
    
    catch_area = data.table(
      ara_id = catchments$ara_id, 
      area = st_area(catchments))
    catch_pops_covs_summed = merge(catch_pops_covs_summed, catch_area, by=c('ara_id'))
    catch_pops_covs_summed[, pop_dens := drop_units(B21BTOT/area)]
    
    
    catchement_population_covariates = catch_pops_covs_summed[, c('ara_id', 'B21BTOT', 'prop_under_20', 'prop_over_65', 'prop_non_ch_eu', 'pop_dens')]
    colnames(catchement_population_covariates) = c('ara_id', 'total_pop', 'prop_under_20', 'prop_over_65', 'prop_non_ch_eu', 'pop_dens')
    
    
    #write.csv(catchement_population_covariates, 'data/catchement_population_covariates.csv')
    
    if(FALSE) { # for now remove employment
    # EMPLOYMENT DATA 
    
    catch_emp_data = st_join(st_transform(emp_data_slim_sf, crs=st_crs(plz_catch_area)), intersections)
    catch_emp_data_no_geom = st_drop_geometry(catch_emp_data )[,c('ara_id', 'B08VZAT')]
    catch_emp_data_dt = data.table(data.frame(catch_emp_data_no_geom))
    catch_emp_data_dt[, B08VZAT := nafill(catch_emp_data_dt$B08VZAT, fill=0)]
    catch_emp_data_dt[, catch_emp_pop := sum(B08VZAT), by=c('ara_id')]
    catch_emp_data_dt_summed = unique(catch_emp_data_dt[, c('ara_id', 'catch_emp_pop')])
    
    
    
    ### Combine 
    
    catch_pop_cov = left_join(catchement_population_covariates, catch_emp_data_dt_summed, by = 'ara_id')
    
    
    }
    
    catch_pop_cov = catchement_population_covariates
    #write.csv(catch_emp_data_dt_summed, 'data/catchement_employment_covariates.csv')
    
    return(catch_pop_cov)
    
  } else if(scale == 'PLZ') {
    
    
    
    plz_emp_data = st_join(st_transform(emp_data_slim_sf, crs=st_crs(plz_area)), plz_area)
    plz_emp_data = data.table(plz_emp_data)
    
    plz_emp_data[, plz_emp_pop := sum(B08VZAT), by=c('PLZ')]
    
    
    select_cols = c('B21BTOT', males, fmale, non_CH_EU, "E_KOORD", "N_KOORD")
    
    pop_data_age = pop_data[, ..select_cols]
    
    pop_data_age_sf = st_as_sf(pop_data_age, coords = c("E_KOORD", "N_KOORD"), 
                               crs = bfs_crs, remove=FALSE)
    
    plz_pop_data_age = data.table(st_join(plz_area, st_transform(pop_data_age_sf, crs=st_crs(plz_area))))
    
    
    plz_age_cols = c('PLZ', 'B21BTOT', males, fmale)
    
    plz_pop_data_age_simp = plz_pop_data_age[, ..plz_age_cols]
    
    plz_pop_data_age_simp_summed = plz_pop_data_age_simp[, lapply(.SD, sum, na.rm=TRUE), by=PLZ, .SDcols=c('B21BTOT', males,fmale) ] 
    
    plz_pop_data_age_simp = unique(plz_pop_data_age_simp_summed)
    
    under_20_cols = grep("^B\\d{2}(BM|BW)(0[1-3])$", names(plz_pop_data_age_simp), value = TRUE)
    over_65_cols = grep("^B\\d{2}(BM|BW)(1[4-9]|2\\d|3[0-9])$", names(plz_pop_data_age_simp), value = TRUE)
    plz_pop_data_age_simp[, prop_under_20 := rowSums(.SD) / B21BTOT, .SDcols = under_20_cols  ]
    plz_pop_data_age_simp[, prop_over_65 := rowSums(.SD) / B21BTOT, .SDcols = over_65_cols  ]
    
    
    plz_eth_cols = c('PLZ','B21BTOT', non_CH_EU)
    
    plz_pop_data_birth = plz_pop_data_age[, ..plz_eth_cols]
    
    plz_pop_data_birth_summed = plz_pop_data_birth[, lapply(.SD, sum, na.rm=TRUE), by=PLZ, .SDcols=c('B21BTOT', non_CH_EU)] 
    
    plz_pop_data_birth = unique(plz_pop_data_birth_summed)
    plz_pop_data_birth[, non_ch_eu := rowSums(.SD), .SDcols = non_CH_EU]
    plz_pop_data_birth[, prop_non_ch_eu := non_ch_eu/B21BTOT]
    plz_pop_data_birth$prop_non_ch_eu = nafill(plz_pop_data_birth$prop_non_ch_eu,fill=0)
    
    plz_population_data = plz_pop_data_age_simp[, c('PLZ', 'B21BTOT', 'prop_under_20', 'prop_over_65')]
    
    plz_pop_cov = merge(merge(plz_population_data, unique(plz_emp_data[,c('PLZ', 'plz_emp_pop')]), on='PLZ'), plz_pop_data_birth[, c('PLZ', 'prop_non_ch_eu')])
    plz_area = rmapshaper::ms_dissolve(plz_area,'PLZ')
    plz_surarea = data.table( PLZ = plz_area$PLZ, 
                              area = st_area(plz_area))
    plz_surarea [, area := sum(area), by = c('PLZ')]
    
    plz_pop_cov = merge(plz_pop_cov, plz_surarea, by=c('PLZ'), how='left')
    plz_pop_cov[, pop_dens:=drop_units(B21BTOT/area)]
    
    #write.csv(drop_na(plz_pop_cov[, -c('area')]), 'data/population_statistics_plz.csv')
    
    return(plz_pop_cov)
    
  } else {
    print('Please set scale to "ARA" or "PLZ"')
  }
  

  
}
