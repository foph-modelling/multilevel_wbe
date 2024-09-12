library(patchwork)


model_names =  c('NUTS2', 'Spatial Cluster', 'Mobility - \nDegree', 'Mobility - \nBetweeness', 'Mobility - \nPage Rank', 'Mobility - \nCommuting Hubs')
names(model_names) = as.character(1:6)


scores_wtbed = readRDS('outputs/Full_Run_wtbed_2/scores.rds')
scores_nobed = readRDS('outputs/Full_Run_nobed/scores_default_lab.rds')


scores_wtbed_all = scores_wtbed$scores[time>40,]
scores_nobed_all =  scores_nobed$scores
scores_wtbed_all[, fit := 'With training period']
scores_nobed_all[, fit := 'No training period']

scores_all = rbind(scores_wtbed_all, scores_nobed_all)
  
  
scores_by_model = summarise_scores(scores_all, by = c("model", "fit"))
scores_by_model_long = melt(scores_by_model, id.vars=c('model', 'fit'), measure.vars = colnames(scores_by_model[,-c('model', 'fit')]), value.name = 'score', variable.name='metric')

scores_by_model_long %>% ggplot() + 
  geom_point(aes(x=score, y=model, color=as.character(model))) + 
  geom_vline(xintercept = 0) + 
  scale_y_continuous(transform = 'reverse', breaks=1:6) + 
  facet_grid(fit~metric,scales = 'free') + 
  scale_color_discrete(name='WWTP selection regime', labels=model_names)+
  theme_minimal()+
  theme(
    legend.position = 'bottom'
  )

ggsave('plots/spatial_analysis/scores_overall.png', width=25, height=9, units='cm')



scores_by_time_wtbed = scores_wtbed$scores_by_time[time>40,]
scores_by_time_wtbed[, time := time - 40 ]
scores_by_time_nobed = scores_nobed$scores_by_time
scores_by_time_wtbed[, fit:='With training period']
scores_by_time_nobed[, fit:='No training period']

scores_by_time_all = rbind(scores_by_time_wtbed, scores_by_time_nobed)


sbt_plot_all = scores_by_time_all %>% ggplot(aes(x=time, y=crps, color=as.character(model))) + 
  geom_point() + 
  ggbump::geom_bump() + 
  facet_wrap(~fit, nrow=2)+
  scale_color_discrete(name='wwwp selection regime', labels=model_names)+
  theme_minimal()



sbt_dens = scores_by_time_all %>% ggplot() + geom_density(aes(y=crps, color = as.character(model), group=model)) + facet_wrap(~fit, nrow=2) +
  scale_color_discrete(guide=F)+
  ylab('')+
  theme_minimal()

sbt_plot_all + sbt_dens + plot_layout(widths= c(4,1), guides='collect')

ggsave('plots/spatial_analysis/scores_by_time.png', width=25, height=15, units='cm')


scores_by_catch = summarise_scores(scores_all, by = c("model", "fit", "ara_id"))
scores_by_catch_long = melt(scores_by_catch, id.vars=c('model', 'fit', "ara_id"), measure.vars = colnames(scores_by_catch[,-c('model', 'fit', 'ara_id')]), value.name = 'score', variable.name='metric')
catchments_scores = merge(catchments, scores_by_catch_long, on=c('ara_id'))


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
  catchments_scores %>% filter(fit=='No training period') %>% filter(metric=='crps') %>% ggplot() + 
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
