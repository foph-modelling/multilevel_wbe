#' ---
#' title: "WBE for SARS-CoV-2 in Switzerland: model development (A3) "
#' author: "Julien Riou"
#' date: "`r Sys.Date()`"
#' params:
#'    controls: controls
#' output:
#'    html_document:
#'      code_folding : hide
#'      number_sections: false
#'      highlight: pygments
#'      theme: cosmo
#'      link-citations: true
#' ---


#+ results="hide", warnings="false", echo="false"
# scp savepoints/savepoint_2023-05-15/controls.rds UBELIX:/storage/homefs/jr18s506/projects/multilevel_wbe/savepoints/savepoint_2023-05-15/.
# scp savepoints/savepoint_2023-05-15/ww1.rds UBELIX:/storage/homefs/jr18s506/projects/multilevel_wbe/savepoints/savepoint_2023-05-15/.
# scp UBELIX:/storage/homefs/jr18s506/projects/multilevel_wbe/savepoints/savepoint_2023-05-15/ma5.3.1.rds savepoints/savepoint_2023-05-15/. 
# scp UBELIX:/storage/homefs/jr18s506/projects/multilevel_wbe/savepoints/savepoint_2023-05-15/ma5.3.2.rds savepoints/savepoint_2023-05-15/. 
# scp UBELIX:/storage/homefs/jr18s506/projects/multilevel_wbe/savepoints/savepoint_2023-05-15/ma5.3.3.rds savepoints/savepoint_2023-05-15/. 
# scp UBELIX:/storage/homefs/jr18s506/projects/multilevel_wbe/savepoints/savepoint_2023-05-15/ma5.4.1.rds savepoints/savepoint_2023-05-15/. 

#source("R/setup.R")

library(data.table)

library(dplyr     )   
library(readr)
library(forcats   )
library(stringr )
library(ggplot2   )
library(tibble  )
library(lubridate )  
library(tidyr)
library(purrr)

library(lubridate)
library(ISOweek)
library(INLA)
library(sf)
library(splines)
library(cowplot)
library(spdep)
library(jsonlite)
library(scales,)
library(units)



source('R/mw_009_spacetime_model_functions.R')
source('R/mw_011_evaluation_functions.R')
source('R/mw_008_load_pop_covars.R')


message("sourced functions")

if( !dir.exists('outputs')){
  dir.create('outputs')
}


save.point = paste0('outputs/last_run_', Sys.time())
dir.create(save.point)
ww1 = readRDS(fs::path(controls$savepoint,"ww1.rds"))
shapes = readRDS(fs::path(controls$savepoint,"shapes.rds"))


message("Loaded data")

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


ww_all = ww_all %>% filter(date > lubridate::ymd(20220210) & date < lubridate::ymd(20220630))
#ww_all = ww_all %>% filter(day1<20)

ww_all = ww_all %>% complete(ara_id, day)

ww_all$day1 = ww_all$day

message("Combined ww data and covariates")
# correspondence table
#corr_all = ww_all %>% 
#  group_by(ara_n,ara_id,ara_name,kt,pop,lab,lab_n,lab2,lab_n2,lab_method,lab_method_n,ara1,ara2,NUTS2_name) %>% 
#  count() %>% 
#  ungroup() 
#corr_all_ara = ww_all %>% 
#  group_by(ara1,ara_name,ara_id,kt,NUTS2_name) %>% 
#  count() %>% 
#  ungroup() 
#saveRDS(corr_all_ara,file=paste0("../",controls$savepoint,"corr_all_ara.rds"))


# Construct spatial model inputs

# convert ara.shapes to North/Easting to correspond to distances in m and extract area based centroids (potential improvement might be to use population weighted centoidds)
catchments = shapes$ara_shp
catchments = st_transform(catchments, 25830)
catchment_centroids = st_centroid(catchments)


message("Loaded shape data")

# merge ww data with centroids to ensure geometries map properly


# extract coordinates - non-sf object 


# Make sure no VL vals are 0 or negative - scale down viral loads to support INLA tractability 
ww_all = ww_all %>% mutate(vl_stand = if_else(vl<=0, 1e-23, vl/mean(na.exclude(ww_all$vl))))
ww_all = data.table(ww_all)[order(day, ara_id), ]

# construct list of coordinates to match the viral load data time series data 
catchment_centroids = merge(catchment_centroids, ww_all, by='ara_id', how='right')
centroid_coords = st_coordinates(catchment_centroids)
colnames(centroid_coords) = c('X1', 'Y1')

message("Prepared inputs")


message("Running INLA model ... (expect a long pause)")

# construct and run inla model in fit_inla_model() - in mw_009_spacetime_model_funtions.r
inla_results = fit_inla_model(wwdata = ww_all, 
                              catchment_centroids = catchment_centroids, 
                              plz_pos = plz_pos,
                              save.point = save.point
                              )


message("Loading spatial data for projection")

plz_pos = get_plz_centroids(crs_required = 25830)
plz_coords = cbind(st_drop_geometry(plz_pos[,c('PLZ')]), st_coordinates(plz_pos) )

plz_area = get_plz_areas(crs_required = 25830)

plz_covariate_matrix = mw_008_load_pop_covars(scale = 'PLZ')
message("Loaded spatial data for projection")

time_steps = seq(1, length(unique(ww_all$day)))

pcoords = data.table()
for(time in time_steps){
  print(dim(pcoords))
  pcoords_date = cbind(plz_coords, time)
  pcoords = rbind(pcoords, pcoords_date)
}
names(pcoords) <- c("PLZ", "x", "y", "time")

message("Generated container for samples")
message("Preparing projection covariates")
pred_coords_covars = merge(unique(pcoords), unique(plz_covariate_matrix), by=c('PLZ'), how='left')
message("Prepared projection covariates")


message("Sampling the INLA model...")

covariates = c('u20', 'o65', 'nec', 'pop_dens')
get_samples_from_inla_model(inla_results = inla_results, 
                            covariates = covariates, 
                            pred_coords_covars = pred_coords_covars, 
                            nsims = 500, 
                            model_dir = save.point)

message(paste0("Sampling complete... outputs saved at ", save.point))

message("Scoring samples...")

scores = score_by_catch(nsims = 500,
                         savepath=save.point,
                         pred_coords_covars = pred_coords_covars, 
                         models = c(''),
                         suffix='', 
                         log_vals=T, 
                         buffer=0)
message("Scores generated")

message("Plotting outputs and saving")

scores$all_catch_res_long[, ':='(pred_mean=mean(prediction), upper=quantile(prediction, 0.95, na.rm=T), lower=quantile(prediction, 0.05, na.rm=T)), by=c('time', 'ara_id', 'model')]

catch_summary = unique(scores$all_catch_res_long[, c('time', 'ara_id', 'model', 'pred_mean', 'upper', 'lower', 'true_value')])

catch_summary %>% ggplot() + 
  
  #geom_rect(aes(xmin = day-0.5, xmax = day+0.5, ymin=-Inf, ymax=Inf, fill=lab ), alpha=0.2)+
  geom_point(aes(x=time, y=true_value), color='black', size=0.2) + 
  geom_line(aes(x=time, y=pred_mean, color=model), alpha=0.7) +
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper, fill=model), alpha=0.2)+
  facet_wrap(~ara_id, ncol=10)+
  coord_cartesian(ylim=c(0,750))+ 
  ylab('Viral load per person') + 
  theme_minimal()#+
  #geom_rect(data = unique(subset(catch_summary %>% select(name, model),(name %in% select_wwtps)  & !(name %in% clust_pop[rankpop==1, ]$ara_id))),# &  model=='all_daily')), 
  #          fill = NA, colour = "limegreen", xmin = -Inf,xmax = Inf,
  #          ymin = -Inf,ymax = Inf, linewidth=2)+
  #geom_rect(data =unique(subset(catch_summary %>% select(name, model),!(name %in% select_wwtps)  & (name %in% clust_pop[rankpop==1, ]$ara_id))),# &  model=='all_daily'))), 
  #          fill = NA, colour = "red", xmin = -Inf,xmax = Inf,
  #          ymin = -Inf,ymax = Inf, linewidth=2)+
  #geom_rect(data = unique(subset(catch_summary %>% select(name, model),(name %in% select_wwtps)  & (name %in% clust_pop[rankpop==1, ]$ara_id))),# &  model=='all_daily')), 
  #          fill = NA, colour = "purple", xmin = -Inf,xmax = Inf,
  #          ymin = -Inf,ymax = Inf, linewidth=2)
#
#scale_y_continuous(trans='log')

ggsave('catchplot.png', height=10, width=15, units='in')

dp_s = copy(ww_all)

index_s <- inla.stack.index(stack = inla_results$stk, tag = "est")$data
dp_s$pred_mean = inla_results$res$summary.fitted.values[index_s, "mean"]
dp_s$pred_lower = inla_results$res$summary.fitted.values[index_s, "0.025quant"]
dp_s$pred_upper = inla_results$res$summary.fitted.values[index_s, "0.975quant"]



ggplot(dp_s[below_lod==0 & below_loq==0, ]) + 
  geom_point(aes(x=day1, y=vl_stand), size=0.3) + 
  geom_line(aes(x=day1, y=pred_mean), color='red', linewidth=0.2, alpha=0.8)+
  geom_ribbon(aes(x=day1, ymin=pred_lower, ymax=pred_upper), fill='red', alpha=0.2)+
  facet_wrap(~ara_id)+
  scale_y_continuous(trans='log')

ggsave('ppp.png', height=10, width=15, units='in')

