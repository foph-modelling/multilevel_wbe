
### save.point = 
ww_all = readRDS(file = paste0(save.point, '/ww_all.rds'))
catchment_centroids = readRDS(paste0(save.point, '/catchment_centroids.rds'))
selction_ids = readRDS('outputs/select_catchments.rds')

selection = 1
ara_ids = selction_ids[[selection]]
ww_tofit = ww_all[!(ara_id %in% ara_ids) & date > lubridate::ymd(20220320),  vl_stand := NA]

# construct list of coordinates to match the viral load data time series data 
catchment_centroids_select = merge(catchment_centroids, ww_tofit, by='ara_id', how='right')
centroid_coords = st_coordinates(catchment_centroids_select)
colnames(centroid_coords) = c('X1', 'Y1')

# construct and run inla model in fit_inla_model() - in mw_009_spacetime_model_funtions.r
inla_results_select = fit_inla_model(wwdata = ww_tofit, 
                                     catchment_centroids = catchment_centroids_select, 
                                     plz_pos = plz_pos,
                                     save.point = save.point, 
)

saveRDS(inla_results_select, paste0(save.point, '/model_', selection, '.rds'))