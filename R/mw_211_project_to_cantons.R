library(data.table)
library(sf)
library(ggplot2)



project_to_cantons = function(savepath, pcoords, start, catchments, models=c('')){
  #Cantonal boundaries
  cantons = st_read('data/spatial/bfs/swissboundaries3d_2023-01_2056_5728.shp/swissBOUNDARIES3D_1_4_TLM_KANTONSGEBIET.shp')
  cantons = st_transform(cantons, 25830)
  
  log_res = TRUE
  covariate_pops_plz = fread('data/population_statistics/population_statistics_plz.csv')
  
  
  
  #generate list of PLZ in each Canton
  plz_pos = get_plz_centroids(crs_required=25830)
  PLZ_cant = st_join(cantons, plz_pos)
  PLZ_cant_list = st_drop_geometry(PLZ_cant)[,c('PLZ', 'NAME')]
  
  
  
  catchments = st_transform(catchments, 25830)
  
  WWVLS_summary = ww_all
  #WWVLS_summary = process_ww_data_daily(WWVLS_summary = data.table(vl_data$long), start = 20220210, end=20220331)
  #WWVLS_summary = merge(WWVLS_summary, ww_all[, c("name", "date", "lab", "method")], by=c('name', 'date'))
  #WWVLS_summary[, vl_mean7d := if_else(vl<=0, 1e-10, vl/1e11)]
  
  
  all_cant_res_long = data.table()
  for (model in models) {
    nsims = 500
    
    res_path = paste0('outputs/samples_new_gam_pcp/', model, '/posterior_predictions_', nsims, '_', model, '.rds')
    inp_path = paste0('outputs/samples_new_gam_pcp/', model, '/model_inputs.rds')
    
    
    res_path = paste0(savepath, '/posterior_predictions_', nsims, '_', model, '.rds')
    inp_path = paste0(savepath, '/model_inputs.rds')
    
    
    sims.pred.t = readRDS(res_path)
    #mod_inps = readRDS(inp_path)
    
    #pcoords = pred_coords_covars
    pcoords = pcoords[order(time, PLZ),]
    
    
    colnames(sims.pred.t) = paste0('sample_', seq(1,nsims))
    
    sim.samples.indexed = pcoords[,c('PLZ', 'time')]
    
    if(log_res==TRUE){
      print('log res true')
      sim.samples.indexed = cbind(sim.samples.indexed, exp(sims.pred.t))
    }else{
      sim.samples.indexed = cbind(sim.samples.indexed, sims.pred.t)
    }
    
    
    
    
    
    
    catchments_meas = merge(catchments, WWVLS_summary, by=c('ara_id'))
    
    
    # load and preprocess the population weight data for cantons 
    plz_cant_populations = data.table(merge(PLZ_cant_list, covariate_pops_plz[, c('PLZ', 'B21BTOT')], by=c('PLZ')))
    
    plz_cant_populations[, weight:= B21BTOT/sum(B21BTOT), by=c('NAME')]
    plz_cant_populations[, weight_test := sum(weight), by = c('NAME')]
    
    plz_cant_populations[, sw_weight := B21BTOT/sum(B21BTOT)]
    
    
    sampcols = paste0('sample_', seq(1,nsims))
    sampcolswt = paste0('samplewt_', seq(1,nsims))
    sampcolscant = paste0('samplecant_', seq(1,nsims))
    
    # merge results with population weights by catchment and calculate weighted VLs 
    weighted_ww_vl_plz_cantons = merge(sim.samples.indexed, plz_cant_populations[,c('PLZ', 'NAME', 'weight')], by=c('PLZ'))
    weighted_ww_vl_plz_cantons[, (sampcolswt) := lapply(.SD, function(x){x*weight}), .SDcols = sampcols]
    weighted_ww_vl_plz_cantons[, (sampcolscant) := lapply(.SD, sum), by = c('NAME', 'time'), .SDcols = sampcolswt]
    
    
    samp_vl_cantons = weighted_ww_vl_plz_cantons[,c('NAME', 'time', ..sampcolscant)]
    
    remove(weighted_ww_vl_plz_cantons)
    samp_vl_cantons_unique = data.table()
    for(n in unique(samp_vl_cantons$NAME)){
      for(t in unique(samp_vl_cantons$t)){
        samp_vl_cantons_unique = rbind(samp_vl_cantons_unique, unique(samp_vl_cantons[NAME==n & time==t,]))
      }
    }
    
    
    samp_true_cantons = samp_vl_cantons_unique
    
    samp_true_cantons_long = melt(samp_true_cantons, id.vars = c('NAME', 'time'), measure.vars = sampcolscant, variable.name = 'sample', value.name = 'prediction')
    
    
    samp_true_cantons_long[, model := model]
    
    all_cant_res_long = rbind(all_cant_res_long, samp_true_cantons_long)
    
  }
  

  reference_date = start
  hosps = fread('../../02_data/bag_covid_19_data_csv_12_September_2023/data/COVID19Hosp_geoRegion.csv')
  
  hosps_ct = hosps[!(geoRegion %in% c('CH', 'CHFL','FL')),]
  geotoname = fread('data/cantoncodes.csv', skip=1)
  colnames(geotoname) = c('ind', 'geoRegion', 'NAME')
  hosps_ct = merge(hosps_ct, geotoname, by=c('geoRegion'))
  
  
  hosps_ct[, day := as.numeric(as.Date(datum) - reference_date) + 1]
  
  
  all_cant_res_long[, ':='(pred_mean=mean(prediction), upper=quantile(prediction, 0.975, na.rm=TRUE), lower=quantile(prediction, 0.025, na.rm=TRUE)), by=c('time', 'NAME', 'model')]
  
  cant_summary = unique(all_cant_res_long[, c('time', 'NAME', 'model', 'pred_mean', 'upper', 'lower')])
  
  cant_summary %>% ggplot() + 
    geom_line(aes(x=time, y=pred_mean, color=as.character(model), group=model), alpha=0.7) +
    geom_ribbon(aes(x=time, ymin=lower, ymax=upper, fill=as.character(model), group=model), alpha=0.2)+
    #geom_line(aes(x=time, y=exp(pred_mean), color=model), alpha=0.7) +
    #geom_ribbon(aes(x=time, ymin=exp(lower), ymax=exp(upper), fill=model), alpha=0.2)+
    geom_point(data=hosps_ct[day %in% cant_summary$time][day %in% cant_summary$time], aes(x=day, y=mean7d*1e6/(4*pop)), color='black', size=0.3)+
    facet_wrap(~NAME, ncol=7)+
    ylab('Viral load per person')+
    coord_cartesian(ylim=c(0,10))+ 
    theme_minimal() + 
    theme(legend.position = 'none')
  ggsave(paste0(save.point, '/cantons_vl_plus_hosps.png'), width = 15, height=10, units='in')
  
  
  
  cases = fread('../../02_data/bag_covid_19_data_csv_12_September_2023/data/COVID19Cases_geoRegion.csv')
  
  cases_ct = cases[!(geoRegion %in% c('CH', 'CHFL','FL')),]
  geotoname = fread('data/cantoncodes.csv', skip=1)
  colnames(geotoname) = c('ind', 'geoRegion', 'NAME')
  cases_ct = merge(cases_ct, geotoname, by=c('geoRegion'))
  
  
  cases_ct[, day := as.numeric(as.Date(datum) - reference_date) + 1]
  
  
 
  cant_summary %>% ggplot() + 
    geom_line(aes(x=time, y=pred_mean, color=as.character(model), group=model), alpha=0.7) +
    geom_ribbon(aes(x=time, ymin=lower, ymax=upper, fill=as.character(model), group=model), alpha=0.2)+
    #geom_line(aes(x=time, y=exp(pred_mean), color=model), alpha=0.7) +
    #geom_ribbon(aes(x=time, ymin=exp(lower), ymax=exp(upper), fill=model), alpha=0.2)+
    geom_point(data=cases_ct[day %in% cant_summary$time], aes(x=day, y=mean7d*8e2/pop), color='black', size=0.3)+
    facet_wrap(~NAME, ncol=7)+
    ylab('Viral load per person')+
    coord_cartesian(ylim=c(0,10))+ 
    theme_minimal() + 
    theme(legend.position = 'none')
  ggsave(paste0(save.point, '/cantons_vl_plus_cases.png'), width = 15, height=10, units='in')
  
  saveRDS(cant_summary, paste0(save.point, '/cant_summary.rds'))
  saveRDS(cases_ct, paste0(save.point, '/cases_ct.rds'))
  saveRDS(hosps_ct, paste0(save.point, '/hosps_ct.rds'))
  saveRDS(all_cant_res_long, paste0(save.point, '/all_cant_res_long.rds'))
  
}





