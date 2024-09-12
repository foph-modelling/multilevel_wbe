pred_coords_covars = readRDS(paste0(save.point, '/pred_coords_covars.rds'))
for(selection in 1:6){
  inla_results = readRDS(paste0(save.point, 'model_', selection, '.rds'))
 
  covariates = c('u20', 'o65', 'nec', 'log_pop_dens')
  get_samples_from_inla_model(inla_results = inla_results, 
                              covariates = covariates, 
                              pred_coords_covars = pred_coords_covars, 
                              nsims = 500, 
                              model_dir = save.point, 
                              model=selection)
}


scores = score_by_catch(nsims = 500,
                        ww_all = ww_all,
                        savepath=save.point,
                        pred_coords_covars = pred_coords_covars, 
                        models = c(0),
                        suffix='', 
                        log_vals=T)


scores$all_catch_res_long[, ':='(pred_mean=mean(prediction), upper=quantile(prediction, 0.95, na.rm=T), lower=quantile(prediction, 0.05, na.rm=T)), by=c('time', 'ara_id', 'model')]

catch_summary = unique(scores$all_catch_res_long[, c('time', 'ara_id', 'model', 'pred_mean', 'upper', 'lower', 'true_value')])

catch_summary %>% ggplot() + 
  
  #geom_rect(aes(xmin = day-0.5, xmax = day+0.5, ymin=-Inf, ymax=Inf, fill=lab ), alpha=0.2)+
  geom_point(aes(x=time, y=true_value), color='black', size=0.2) + 
  geom_line(aes(x=time, y=pred_mean, color=as.character(model), group=model), alpha=0.7) +
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper, fill=as.character(model), group=model), alpha=0.2)+
  facet_wrap(~ara_id)+
  scale_color_discrete(name='Model', labels=model_names)+ 
  scale_fill_discrete(name='Model', labels=model_names)+ 
  coord_cartesian(ylim=c(0,5))+ 
  ylab('Viral load per person') + 
  theme_minimal()#+


scores_by_model = summarise_scores(scores$scores[time>40,], by = c("model"))
scores_by_model_long = melt(scores_by_model, id.vars=c('model'), measure.vars = colnames(scores_by_model[,-c('model')]), value.name = 'score', variable.name='metric')

scores_by_model_long %>% ggplot() + 
  geom_point(aes(x=score, y=model, color=as.character(model))) + 
  geom_vline(xintercept = 0) + 
  scale_y_continuous(transform = 'reverse', breaks=1:6) + 
  facet_wrap(~metric, nrow=1, scales = 'free') + 
  theme_minimal()



scores$scores_by_time %>% ggplot() + geom_line(aes(x=time, y=crps, color=as.character(model)))


scores$scores_by_catch %>% 
  ggplot(aes(x=model, y=crps, color=ara_id)) + 
  geom_point()+
  ggbump::geom_bump() + 
  facet_wrap(~ara_id, ncol=5)+
  theme(legend.position='none')


sbt_plot = scores$scores_by_time %>% ggplot(aes(x=time, y=crps, color=as.character(model))) + 
  geom_point() + 
  ggbump::geom_bump() + 
  theme_bw()


scores$scores_by_time %>% ggplot() + geom_density(aes(x=mad, color = as.character(model), group=model))

saveRDS(scores, paste0(save.point, '/scores.rds'))


project_to_cantons(save.point, pcoords = pred_coords_covars, start = min(ww_all$date, na.rm = T), catchments = catchments, models=c(1, 2, 3, 4, 5, 6) )
