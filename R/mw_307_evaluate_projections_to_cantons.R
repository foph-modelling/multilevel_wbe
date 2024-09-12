cant_summary = readRDS('outputs/Full_Run_nobed_all_nolab/cant_summary.rds')
cant_summary_all = readRDS('outputs/Full_Run_all/cant_summary.rds')

cases_ct = readRDS(paste0(save.point, '/cases_ct.rds'))
hosps_ct = readRDS(paste0(save.point, '/hosps_ct.rds'))

cases_ct = cases_ct[day %in% 1:78,]
hosps_ct = hosps_ct[day %in% 1:78,]

x = unlist(cases_ct[NAME=='Genève',]$mean7d)
y = unlist(cant_summary[NAME=='Genève' & model==5,]$pred_mean)

trans_ent = data.table()
for(name in unique(cases_ct$NAME)){
  for(mod in 1:6){
    x = unlist(cases_ct[NAME==name,]$mean7d)
    y = unlist(cant_summary[NAME==name & model==mod,]$pred_mean)
    plot((y-sd(y))/mean(y)) + lines((x-sd(x))/mean(x))
    trans_ent = rbind(trans_ent,
      data.table(NAME = name, 
                 model = mod, 
                 TE = RTransferEntropy::calc_te(y[20:78], x[20:78], lx = 5, ly = 5),
                 SR = DescTools::SpearmanRho(x,y), 
                 DTW = dtwclust::dtw_basic(as.matrix((x - sd(x))/mean(x)) , as.matrix((y - sd(y))/mean(y))),
                 ccf = max(ccf((x - sd(x))/mean(x) , (y - sd(y))/mean(y))[[1]]), 
                 mad = mean(scoringutils::mad_sample(matrix(all_cant_res_long[NAME==name & model==mod,]$prediction, nrow=78)))))
    
  }
}


trans_ent_all = data.table()
for(name in unique(cases_ct$NAME)){
    x = unlist(cases_ct[NAME==name,]$mean7d)
    y = unlist(cant_summary_all[NAME==name,]$pred_mean)
    trans_ent_all = rbind(trans_ent_all,
                      data.table(NAME = name, 
                                 model = mod, 
                                 TE = RTransferEntropy::calc_te(y[20:78], x[20:78], lx = 5, ly = 5),
                                 SR = DescTools::SpearmanRho(x,y), 
                                 DTW = dtwclust::dtw_basic(as.matrix((x - sd(x))/mean(x)) , as.matrix((y - sd(y))/mean(y))),
                                 ccf = max(ccf((x - sd(x))/mean(x) , (y - sd(y))/mean(y))[[1]]), 
                                 mad = mean(scoringutils::mad_sample(matrix(all_cant_res_long[NAME==name & model==mod,]$prediction, nrow=78)))))
    
  }





trans_ent_hosp = data.table()
for(name in unique(cases_ct$NAME)){
  for(mod in 1:6){
    x = unlist(hosps_ct[NAME==name,]$mean7d)
    y = unlist(cant_summary[NAME==name & model==mod,]$pred_mean)
    trans_ent_hosp = rbind(trans_ent_hosp,
                      data.table(NAME = name, 
                                 model = mod, 
                                 TE = RTransferEntropy::calc_te(y[20:78], x[20:78], lx = 5, ly = 5),
                                 SR = DescTools::SpearmanRho(x,y), 
                                 DTW = dtwclust::dtw_basic(as.matrix((x - sd(x))/mean(x)) , as.matrix((y - sd(y))/mean(y))), 
                                 ccf = max(ccf((x - sd(x))/mean(x) , (y - sd(y))/mean(y))[[1]]), 
                                 mad = mean(scoringutils::mad_sample(matrix(all_cant_res_long[NAME==name & model==mod,]$prediction, nrow=78)))))
                      
    
  }
}


trans_ent_hosp_all = data.table()
for(name in unique(cases_ct$NAME)){
    x = unlist(hosps_ct[NAME==name,]$mean7d)
    y = unlist(cant_summary_all[NAME==name,]$pred_mean)
    trans_ent_hosp_all = rbind(trans_ent_hosp_all,
                           data.table(NAME = name, 
                                      model = mod, 
                                      TE = RTransferEntropy::calc_te(y[20:78], x[20:78], lx = 5, ly = 5),
                                      SR = DescTools::SpearmanRho(x,y), 
                                      DTW = dtwclust::dtw_basic(as.matrix((x - sd(x))/mean(x)) , as.matrix((y - sd(y))/mean(y))), 
                                      ccf = max(ccf((x - sd(x))/mean(x) , (y - sd(y))/mean(y))[[1]]), 
                                      mad = mean(scoringutils::mad_sample(matrix(all_cant_res_long[NAME==name & model==mod,]$prediction, nrow=78)))))
    
    
  }




trans_ent[,mean_mad:= mean(mad), by = c('model')]
trans_ent[,mean_DTW:= mean(DTW), by = c('model')]
trans_ent[,mean_ccf:= mean(ccf), by = c('model')]
trans_ent[,mean_SR:= mean(SR), by = c('model')]


trans_ent_hosp[,mean_mad:= mean(mad), by = c('model')]
trans_ent_hosp[,mean_DTW:= mean(DTW), by = c('model')]
trans_ent_hosp[,mean_ccf:= mean(ccf), by = c('model')]


trans_ent_all[,mean_mad:= mean(mad), by = c('model')]
trans_ent_all[,mean_DTW:= mean(DTW), by = c('model')]
trans_ent_all[,mean_ccf:= mean(ccf), by = c('model')]
trans_ent_all[,mean_SR:= mean(SR), by = c('model')]


trans_ent_hosp_all[,mean_mad:= mean(mad), by = c('model')]
trans_ent_hosp_all[,mean_DTW:= mean(DTW), by = c('model')]
trans_ent_hosp_all[,mean_ccf:= mean(ccf), by = c('model')]


ggplot() + 
  geom_point(data=trans_ent, aes(y=DTW, x=mad, color=as.character(NAME)))+ 
  geom_point(data = unique(trans_ent[,c('model', 'mean_mad', 'mean_DTW')]), aes(x=mean_mad, y=mean_DTW))+
  geom_point(data=unique(trans_ent_all[,c('model', 'mean_mad', 'mean_DTW')][,-c('model')]), aes(y=mean_DTW, x=mean_mad), color='red', shape=3, size=5)+ 
  geom_vline(xintercept=0.5, linetype=2) + 
  geom_hline(yintercept=35, linetype=2) + 
  xlim(c(0,1))+
  scale_color_discrete(name='Canton')+
  facet_wrap(~model, nrow=3, labeller = labeller(model=model_names))+
  theme_minimal()

ggsave('plots/spatial_analysis/score_against_cases.png', width=20, height=15, units='cm')



ggplot() + 
  geom_point(data=trans_ent_hosp, aes(y=DTW, x=mad, color=as.character(NAME)))+ 
  geom_point(data = unique(trans_ent_hosp[,c('model', 'mean_mad', 'mean_DTW')]), aes(x=mean_mad, y=mean_DTW))+
  geom_point(data=unique(trans_ent_hosp_all[,c('model', 'mean_mad', 'mean_DTW')][,-c('model')]), aes(y=mean_DTW, x=mean_mad), color='red', shape=3, size=5)+ 
  geom_vline(xintercept=0.5, linetype=2) + 
  geom_hline(yintercept=40, linetype=2) + 
  xlim(c(0,1))+
  facet_wrap(~model, nrow=3, labeller = labeller(model=model_names))+
  scale_color_discrete(name='Canton')+
  theme_minimal()

ggsave('plots/spatial_analysis/score_against_hosps.png', width=20, height=15, units='cm')

