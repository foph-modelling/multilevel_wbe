#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load population from bfs hectare data and aggregate to plz and ARA
# creation: jmunday
# init date: 2023-08-16
#:::::::::::::::::::::::::::::


mw_008b_load_pop_covars = function(scale='ARA') {
  
  sf_use_s2(FALSE)
  swissboundaries_BFS_NATION = st_read('data/spatial/bfs/swissboundaries3d_2023-01_2056_5728.shp/swissBOUNDARIES3D_1_4_TLM_LANDESGEBIET.shp')
  bfs_crs = st_crs(swissboundaries_BFS_NATION)
  
  # catchments
  catchments = read_sf('data/spatial/230214_ARA_BAG/230214_ARA_BAG.shp')
  catchments = st_make_valid(catchments)
  catchment_crs = st_crs(catchments)
  if(FALSE) {
    ggplot(catchments) + geom_sf()
  }
  
  # population by hectare
  pop_data = data.table::fread('data/population_statistics/STATPOP2021.csv')
  males = names(pop_data)[grepl('^B\\d{2}(BM)(0|1)', names(pop_data))] 
  fmale = names(pop_data)[grepl('^B\\d{2}(BW)(0|1)', names(pop_data))] 
  non_CH_EU = c('B21B28', 'B21B29', 'B21B30')
  select_cols = c('B21BTOT', males, fmale, non_CH_EU, "E_KOORD", "N_KOORD")
  pop_data_slim = pop_data[,..select_cols]
  pop_data_slim_sf = st_as_sf(pop_data_slim, coords = c("E_KOORD", "N_KOORD"), 
                              crs = bfs_crs, remove=FALSE)
  if(FALSE) {
    ggplot(pop_data_slim_sf) + geom_sf(aes(colour=B21BTOT),alpha=.3,size=.2) + scale_colour_viridis_c()
  }
  
  ## Load employment data (Total FTE employment by hectare (2020) (100mx100m squares))
  emp_data = data.table::fread('data/employment_statistics/STATENT_2020.csv',fill = TRUE)
  emp_data_slim = emp_data[,c('N_KOORD', 'E_KOORD', 'B08VZAT')]
  emp_data_slim_sf = st_as_sf(emp_data_slim, coords = c("E_KOORD", "N_KOORD"), 
                              crs = bfs_crs, remove=FALSE)
  if(FALSE) {
    ggplot(emp_data_slim_sf) + geom_sf(aes(colour=B08VZAT),alpha=.3,size=.2) + scale_colour_viridis_c()
  }
  
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
    
  
  
}
