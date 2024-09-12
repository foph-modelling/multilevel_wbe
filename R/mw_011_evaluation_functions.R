library(scoringutils)
library(data.table)



score_by_catch = function(nsims = 500,
                          savepath,
                          ww_all,
                          start='',
                          pred_coords_covars,
                          covariate_matrix,
                          models = c(''),
                          suffix='', 
                          log_vals=F, 
                          buffer=0, 
                          use_default_lab_effect = F){
  if(suffix != ''){
    suffix=paste0('_', suffix)
  }
  
  all_catch_res_long = data.table()
  for (model in models) {
    
    res_path = paste0(savepath, '/', start, 'posterior_predictions_', nsims, '_', model, suffix, '.rds')
    inp_path = paste0(savepath, '/', start, 'model_inputs.rds')
    lab_path = paste0(savepath, '/', start, 'lab_effect_', nsims, '_', model, suffix, '.rds')
    
    if(use_default_lab_effect){
      lab_path = paste0(savepath, '/', start, 'lab_effect_', nsims, '_default.rds')
    }
    
    message(paste0('loading data from ', res_path))
    
    sims.pred.t = readRDS(res_path)
    #mod_inps = readRDS(inp_path)
    lab_effect = readRDS(lab_path)
    
    pcoords = pred_coords_covars
    
    pcoords = pcoords[,c("PLZ", "x", "y", "time")]
    
    print(dim(pcoords)[1]/length(time_steps))
    
    colnames(sims.pred.t) = paste0('sample_', seq(1,nsims))
    
    sim.samples.indexed = pcoords[,c('time', 'PLZ')]
    
    sim.samples.indexed = sim.samples.indexed[order(time, PLZ), ]
    sim.samples.indexed = cbind(sim.samples.indexed, sims.pred.t)
    
    
    catchments_meas = merge(catchments, ww_all, by=c('ara_id'))
    
    
    plz_pop_weights_catchments = fread('data/population_statistics/population_by_plz_in_catchment.csv')
    
    
    
    plz_pop_weights_catchments = plz_pop_weights_catchments %>% rename(ara_id = id)
    
    plz_pop_weights_catchments[, weight := plz_catch_pop / sum(plz_catch_pop), by = c('ara_id')]
    
    plz_pop_weights_catchments[, weight_test := sum(weight), by = c('ara_id')]
    
    sampcols = paste0('sample_', seq(1,nsims))
    sampcolswt = paste0('samplewt_', seq(1,nsims))
    sampcolscatch = paste0('samplecatch_', seq(1,nsims))
    
    # merge results with population weights by catchment and calculate weighted VLs 
    weighted_ww_vl_plz_catchments = merge(sim.samples.indexed, plz_pop_weights_catchments[,c('PLZ', 'ara_id', 'weight')], by=c('PLZ'))
    weighted_ww_vl_plz_catchments[, (sampcolswt) := lapply(.SD, function(x){x*weight}), .SDcols = sampcols]
    weighted_ww_vl_plz_catchments[, (sampcolscatch) := lapply(.SD, sum), by = c('ara_id', 'time'), .SDcols = sampcolswt]
    
    
    samp_vl_catchments = weighted_ww_vl_plz_catchments[,c('ara_id', 'time', ..sampcolscatch)]
    
    
    
    
    
    remove(weighted_ww_vl_plz_catchments)
    samp_vl_catchment_unique = data.table()
    for(n in unique(samp_vl_catchments$ara_id)){
      print(n)
      for(t in unique(samp_vl_catchments$t)){
        samp_vl_catchments_ara = unique(samp_vl_catchments[ara_id==n & time==t,])
        if(length(lab_effect[lab_method_ara[ara_id == n,]$lab_method,]) != 0){
          
          samp_vl_catchments_ara = data.table(t(c(ara_id = n, time = t, unlist(samp_vl_catchments_ara[,-c('ara_id', 'time')] + lab_effect[lab_method_ara[ara_id == n,]$lab_method,]))))
        }
        samp_vl_catchment_unique = rbind(samp_vl_catchment_unique, samp_vl_catchments_ara)
       
      }
    }
    
    samp_vl_catchment_unique[, ara_id := as.character(ara_id)]
    true_values = samp_vl_catchment_unique[,c('ara_id', 'time')]
    
    
    #ww_all[vl_stand > quantile(ww_all$vl_stand, 0.975, na.rm=T), vl_stand := quantile(ww_all$vl_stand, 0.975, na.rm=T)]

    true_values = merge(true_values, ww_all[, c('ara_id', 'day', 'vl_smooth')] %>% mutate(time=day - min(ww_all$day)), on=c('ara_id', 'time'), all.y=T)
    
    
    #mask = which(true_values$name == n & !is.na(true_values$vl_mean7d_smooth))
    
    #scoringutils::crps_sample(true_values[mask,]$vl_mean7d_smooth, as.matrix(samp_vl_catchment_unique[,-c('name', 'time')][mask,]) )
    
    
    samp_true_catchment = merge(true_values, samp_vl_catchment_unique, on=c('time', 'ara_id'))
    
    samp_true_catchment_long = melt(samp_true_catchment, id.vars = c('ara_id', 'time', 'vl_smooth'), measure.vars = sampcolscatch, variable.name = 'sample', value.name = 'prediction')
    
    #samp_true_catchment_long[, vl_smooth := zoo::rollmean(vl_stand, 7, na.rm=T, fill = T), by=c('ara_id')]
    
    samp_true_catchment_long[, true_value:=vl_smooth]
    
    #if(log_vals==T){
    #  samp_true_catchment_long[, true_value := log(true_value + 1e-23)]
    #}
    
    if(log_vals==T){
      samp_true_catchment_long[, prediction := exp(prediction)]
    }
    
    
    samp_true_catchment_long[, model := model]
    
    all_catch_res_long = rbind(all_catch_res_long, samp_true_catchment_long[,-c('vl_stand', 'vl_smooth')])
    
  }
  
  
  scores = scoringutils::score(all_catch_res_long[time > buffer, ])
  
  scores_by_catch = scores %>% scoringutils::summarise_scores(by = c('ara_id', 'model'))
  scores_by_time = scores %>% scoringutils::summarise_scores(by = c('time', 'model'))
  
  
  list(scores=scores, scores_by_catch=scores_by_catch, scores_by_time=scores_by_time, all_catch_res_long= all_catch_res_long)
  
}
