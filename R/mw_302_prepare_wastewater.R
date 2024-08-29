source("R/setup.R")



if( !dir.exists('outputs')){
  dir.create('outputs')
}


save.point = paste0('outputs/last_run_select', Sys.time())
dir.create(save.point)
ww1 = readRDS(fs::path(controls$savepoint,"ww1.rds"))
shapes = readRDS(fs::path(controls$savepoint,"shapes.rds"))

#' 
#' We now model all ARAs together. 
#' 
#' # Switzerland
#' 
#' 
#' 
#' 
ww_all = ww1 %>%
  # log
  dplyr::mutate(logvl=log(vl)) %>%
  mutate(vl=if_else(vl==0, 1, vl)) %>%
  # create indexes for INLA
  dplyr::mutate(day1=day,
                ara1=as.numeric(as.factor(ara_n)),
                ara2=ara1) %>% 
  # group KLBS and KLZH (otherwise not identifiable)
  dplyr::mutate(lab2=if_else(lab=="KLBS","KLZH",lab),
                lab_n2=as.numeric(as.factor(lab2))) %>% 
  # group lab and method
  dplyr::mutate(lab_method=paste0(lab2,"_",method),
                lab_method_n=as.numeric(as.factor(lab_method)))
saveRDS(ww_all,file=paste0(save.point,"/ww_all.rds"))


ww_all = ww_all %>% filter(date > lubridate::ymd(20220210) & date < lubridate::ymd(20220430))
#ww_all = ww_all %>% filter(day1<20)

ww_all = ww_all %>% complete(ara_id, day)

ww_all$day1 = ww_all$day

ww_all$log_pop_dens = log(ww_all$pop_dens)


# correspondence table
corr_all = ww_all %>% 
  group_by(ara_n,ara_id,ara_name,kt,pop,lab,lab_n,lab2,lab_n2,lab_method,lab_method_n,ara1,ara2,NUTS2_name) %>% 
  count() %>% 
  ungroup() 
corr_all_ara = ww_all %>% 
  group_by(ara1,ara_name,ara_id,kt,NUTS2_name) %>% 
  count() %>% 
  ungroup() 
#saveRDS(corr_all_ara,file=paste0("../",controls$savepoint,"corr_all_ara.rds"))

if(!controls$rerun_models) {
  mw_100_desc_table(ww_all) %>%
    dplyr::mutate(across(everything(),as.character)) %>%
    tidyr::gather() %>%
    flextable::flextable(cwidth=c(4,4))
}

# Construct spatial model inputs

# convert ara.shapes to North/Easting to correspond to distances in m and extract area based centroids (potential improvement might be to use population weighted centoidds)
catchments = shapes$ara_shp
catchments = st_transform(catchments, 25830)
catchment_centroids = st_centroid(catchments)

# merge ww data with centroids to ensure geometries map properly


# extract coordinates - non-sf object 


# Make sure no VL vals are 0 or negative - scale down viral loads to support INLA tractability 
ww_all = ww_all %>% mutate(vl_stand = if_else(vl<=0, 1e-23, vl/mean(na.exclude(ww_all$vl))))
ww_all = data.table(ww_all)[order(day, ara_id),]

saveRDS(ww_all, file = paste0(save.point, '/ww_all.rds'))
saveRDS(catchment_centroids, file = paste0(save.point, '/catchment_centroids.rds'))



plz_pos = get_plz_centroids(crs_required = 25830)
plz_coords = cbind(st_drop_geometry(plz_pos[,c('PLZ')]), st_coordinates(plz_pos) )

plz_area = get_plz_areas(crs_required = 25830)

plz_covariate_matrix = mw_008_load_pop_covars(scale = 'PLZ')


time_steps = seq(1, length(unique(ww_all$day)))

pcoords = data.table()
for(time in time_steps){
  print(dim(pcoords))
  pcoords_date = cbind(plz_coords, time)
  pcoords = rbind(pcoords, pcoords_date)
}
names(pcoords) <- c("PLZ", "x", "y", "time")
pred_coords_covars = merge(unique(pcoords), unique(plz_covariate_matrix), by=c('PLZ'), how='left')


saveRDS(pred_coords_covars, paste0(save.point, '/pred_coords_covars.rds'))

