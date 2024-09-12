library(patchwork)


model_names =  c('NUTS2', 'Spatial Cluster', 'Mobility - \nDegree', 'Mobility - \nBetweeness', 'Mobility - \nPage Rank', 'Mobility - \nCommuting Hubs')
names(model_names) = as.character(1:6)


scores_all = readRDS('outputs/Full_Run_all/scores.rds')
scores_sl = readRDS('outputs/Full_Run_nobed_all_nolab/scores.rds')


scores_all_all = scores_all$scores
scores_sl_all =  scores_sl$scores
scores_all_all[, fit := 'all']
scores_sl_all[, fit := 'Select']

scores_both = rbind(scores_all_all, scores_sl_all)
  
  
scores_by_model = summarise_scores(scores_both, by = c("model", "fit"))
scores_by_model_long = melt(scores_by_model, id.vars=c('model', 'fit'), measure.vars = colnames(scores_by_model[,-c('model', 'fit')]), value.name = 'score', variable.name='metric')

metrics_plot = c('mad', 'bias', 'crps', 'ae_median')
 
ggplot() + 
  geom_point(data = scores_by_model_long[metric %in% metrics_plot,], aes(x=score, y=model, color=as.character(model))) + 
  geom_vline(scores_by_model_long[fit=='all' & metric %in% metrics_plot,], mapping = aes(xintercept = score), color='red') + 
  geom_vline(xintercept = 0, color='black')+
  scale_y_continuous(transform = 'reverse', breaks=1:6) + 
  facet_grid(~metric,scales = 'free') + 
  scale_color_discrete(name='WWTP selection regime', labels=c('All WWTPs', as.vector(model_names)))+
  theme_minimal()+
  theme(
    legend.position = 'bottom'
  )

ggsave('plots/spatial_analysis/scores_overall.png', width=15, height=9, units='cm')



scores_by_time_all = data.table(scores_all$scores_by_time)

scores_by_time_sl = scores_sl$scores_by_time
scores_by_time_all[, fit:='all']
scores_by_time_sl[, fit:='Select']

scores_by_time_all = rbind(scores_by_time_all, scores_by_time_sl)


sbt_plot_all = scores_by_time_all %>% ggplot(aes(x=time, y=crps, color=as.character(model))) + 
  geom_point() + 
  ggbump::geom_bump() + 
  #facet_wrap(~fit, nrow=2)+
  scale_color_discrete(name='wwwp selection regime', labels=c('All WWTPs', as.vector(model_names)))+
  theme_minimal()



sbt_dens = scores_by_time_all %>% ggplot() + geom_density(aes(y=crps, color = as.character(model), group=model))  +
  scale_color_discrete(guide=F)+
  ylab('')+
  xlab('')+
  theme_minimal()

sbt_plot_all + sbt_dens + plot_layout(widths= c(4,1), guides='collect')

ggsave('plots/spatial_analysis/scores_by_time.png', width=20, height=9, units='cm')


scores_by_catch = summarise_scores(scores_both, by = c("model", "fit", "ara_id"))
scores_by_catch_long = melt(scores_by_catch, id.vars=c('model', 'fit', "ara_id"), measure.vars = colnames(scores_by_catch[,-c('model', 'fit', 'ara_id')]), value.name = 'score', variable.name='metric')
catchments_scores = merge(catchments, scores_by_catch_long, on=c('ara_id'))

select_shapes = shapes$ara_shp %>% mutate(selection = 0) %>% filter(ara_id == 0)
for(name in 1:length(selction_ids)){
  selection_n = shapes$ara_shp %>% filter(ara_id %in% selction_ids[[name]]) %>% mutate(selection=name)
  select_shapes = rbind(select_shapes, selection_n)
}

score_maps = 
catchments_scores  %>% filter(metric=='crps') %>% ggplot() + 
  geom_sf(data=shapes$canton_shp, fill='lightgray', color='white', alpha=0.5)+
  geom_sf(aes(fill=score), color=NA) + 
  geom_sf(data=select_shapes, color='green') + 
  geom_sf(data=shapes$see_shp, fill='lightblue', color=NA)+
  facet_grid(model~fit, labeller = labeller(model=model_names)) + 
  scale_fill_viridis_c(name="CRPS") + 
  theme_minimal() + 
  theme(legend.position = "bottom", 
        panel.grid=element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
score_maps = 
  catchments_scores %>% filter(fit=='Select') %>% filter(metric=='crps') %>% ggplot() + 
  geom_sf(data=shapes$canton_shp, fill='lightgray', color='white', alpha=0.5)+
  geom_sf(aes(fill=score), color=NA) + 
  geom_sf(data=select_shapes, color='red', fill=NA) + 
  geom_sf(data=shapes$see_shp, fill='lightblue', color=NA)+
  facet_wrap(~model, labeller = labeller(model=model_names)) + 
  scale_fill_viridis_c(name="CRPS") + 
  theme_minimal() + 
  theme(legend.position = "bottom", 
        panel.grid=element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

ggsave('plots/spatial_analysis/score_maps.png', width=20, height=15, units='cm')




