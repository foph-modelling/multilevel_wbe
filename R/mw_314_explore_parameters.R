

fixed_posts = inla_results$res$summary.fixed
hyper_par = inla_results$res$summary.hyperpar
lab_method = inla_results$res$summary.random$lab_method

fixed_posts = fixed_posts %>% mutate(ID=rownames(fixed_posts))
hyper_par = hyper_par %>% mutate(ID:=rownames(hyper_par))

covars_dt = rbind(fixed_posts %>% select(-'kld'), hyper_par)

lab_method %>% ggplot() + 
  geom_pointrange(aes(x=mean, xmin=`0.025quant`, xmax=`0.975quant`, y=ID, color=ID)) + 
  geom_vline(xintercept = 0) + 
  ylab('lab') + 
  xlab('effect (log)') + 
  scale_color_discrete(name='')+
  theme_minimal()

fixed_posts %>% filter(ID != 'b0') %>% ggplot() + 
  geom_pointrange(aes(x=mean, xmin=`0.025quant`, xmax=`0.975quant`, y=ID, color=ID)) + 
  geom_vline(xintercept = 0) + 
  ylab('covariate') + 
  xlab('effect (log)') + 
  scale_color_discrete(name='')+
  theme_minimal()

hyper_par %>% filter(ID == 'Range for s') %>% ggplot() + 
  geom_pointrange(aes(x=mean, xmin=`0.025quant`, xmax=`0.975quant`, y=ID, color=ID)) + 
  geom_vline(xintercept = 0) + 
  ylab('') + 
  xlab('meters') + 
  scale_color_discrete(name='')+
  theme_minimal()
