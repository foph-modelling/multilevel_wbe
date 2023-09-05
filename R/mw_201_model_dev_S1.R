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
if(!exists("controls")) controls = readRDS(fs::path("../savepoints/savepoint_2023-05-15/controls.rds"))
source("setup.R")
ww1 = readRDS(fs::path("../",controls$savepoint,"ww1.rds"))
shapes = readRDS(fs::path("../",controls$savepoint,"shapes.rds"))

#' 
#' We now model all ARAs together. 
#' 
#' # Switzerland
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
saveRDS(ww_all,file=paste0("../",controls$savepoint,"ww_all.rds"))

ww_all = ww_all %>% filter(day1<60)

ww_all = ww_all %>% complete(ara_id, day)

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
catchment_centroids = merge(catchment_centroids, ww_all, by='ara_id', how='right')

# extract coordinates - non-sf object 
centroid_coords = st_coordinates(catchment_centroids)

colnames(centroid_coords) = c('X1', 'Y1')

# merge the coordinates back into the ww data 
ww_all = cbind(data.table(centroid_coords), ww_all)


ww_all = ww_all %>% mutate(vl_stand=vl/1e11)

ww_all = ww_all[order(day, ara_id), ]

# construct mesh for the SPDE
max.edge = diff(range(ww_all$Y1))/(3*5)
mesh <- inla.mesh.2d(
  loc = ww_all[,c('X1', 'Y1')], max.edge = max.edge * c(1,2), cutoff=1000
)
plot(mesh)
points(centroid_coords, col = "red")


# generate stochastic pdes for inla model
spde <- inla.spde2.pcmatern(
  mesh = mesh, alpha = 2,
  prior.range = c(10000, 0.01), # P(range < 10000) = 0.01
  prior.sigma = c(3, 0.01) # P(sigma > 3) = 0.01
)


# construct spatial component of model as SPDE 


#group <- 1:length(ww_all$day)
#timesn <- length(unique(group))
indexs <- inla.spde.make.index("s",
                               n.spde = spde$n.spde)


A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(ww_all[,c('X1', 'Y1')]))

# put data into INLA stack for convinience



stk.e = 
  inla.stack(
  tag = "est",
  data = list(y = ww_all$vl),
  A = list(1,A),
  effects = list(data.frame(below_loq=ww_all$below_loq, 
                            below_lod=ww_all$below_lod, 
                            lab_method=ww_all$lab_method,
                            day=ww_all$day, 
                            weekend=ww_all$weekend,
                            ara1=ww_all$ara1, 
                            day1=ww_all$day1, 
                            ara2=ww_all$ara2), s = indexs)
)




#rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.9)))


# model formula 

formula = y ~ -1 +
  f(below_loq,model="iid") +
  f(below_lod,model="iid") +
  f(day,model="rw2", scale.model=TRUE, constr=TRUE,
    hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
  f(hol,model="linear",mean.linear=0,prec.linear=.2) +
  f(ara1,model="iid") +
  f(day1,model="rw1", scale.model=TRUE, constr=TRUE,
    group=ara2, control.group=list(model="iid"),
    hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) + 
  f(s, model=spde)



#' ## Model A5.3.1: unique national trend
#' 
#' We directly apply model A5.3 to all ARAs. We assume one temporal trend for Switzerland, and each ARA is free to deviate from it independently.
#' 
#+ ma5.3.1, fig.width=8, fig.height=12,  R.options = list(width = 1000)
if(controls$rerun_models) {
  ms5.3.1 = inla(formula,
                 family='gamma',
                 data = inla.stack.data(stk.e),
                 control.compute = list(config=TRUE),
                 control.inla=list(cmin=0),
                 control.predictor = list(
                   compute = TRUE,
                   A = inla.stack.A(stk.e)
                 ),
                 safe = TRUE
  )
  saveRDS(ms5.3.1,file=paste0("../",controls$savepoint,"ma5.3.1.rds"))
} else {
  ms5.3.1 = readRDS(file=paste0("../",controls$savepoint,"ma5.3.1.rds"))
}
summary(ma5.3.1)
summary_exp_vl(ma5.3.1,pars="lab|method|hol|weekend")
ppp_vl_ara(ww_all,ma5.3.1)
#+ ma5.3.1b, fig.width=8, fig.height=4,  R.options = list(width = 1000)
avg_time_trend(ww_all,ma5.3.1)
#+ ma5.3.1c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_130_map_relative_vl(ma5.3.1,corr_all_ara,shapes)
#+ ma5.3.1d, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_131_map_deviation_from_average(ma5.3.1,corr_all_ara,ww_all,shapes,12)



dp_s = copy(ww_all)

index_s <- inla.stack.index(stack = stk.e, tag = "est")$data
dp_s$pred_mean = ms5.3.1$summary.fitted.values[index_s, "mean"]
dp_s$pred_lower = ms5.3.1$summary.fitted.values[index_s, "0.025quant"]
dp_s$pred_upper = ms5.3.1$summary.fitted.values[index_s, "0.975quant"]


index_e <- inla.stack.index(stack = stk.e, tag = "est")$effects


ma5.3.1$summary.fitted.values[index_e,]


ww_all = cbind(ww_all, ma5.3.1$summary.fitted.values)

ggplot(dp_s[below_lod==0 & below_loq==0, ]) + 
  geom_point(aes(x=day1, y=vl), size=0.3) + 
  geom_line(aes(x=day1, y=pred_mean), color='red', linewidth=0.2, alpha=0.8)+
  geom_ribbon(aes(x=day1, ymin=pred_lower, ymax=pred_upper), fill='red', alpha=0.2)+
  facet_wrap(~ara_id)+
  scale_y_continuous(trans='log')


#' 
#' ## Model A5.3.2: effect of lab and method change
#' 
#' We add a covariate to measure the effect of the lab and of changes in methodology. To allow identifiability, we have to make sure that there are multiple ARAs per laboratory, so we group together KLBS (only 1 ARA) and KLZH (12 ARAs). The reference lab is ALTGR.
#'  
#+ ma5.3.2a, fig.width=8, fig.height=12,  R.options = list(width = 1000)
#+ 
#+ 
#+ 
#+ 

ww_all = ww_all  %>% arrange(ara_id, day )

group <- ww_all$day+1
timesn <- length(unique(group))
indexs <- inla.spde.make.index("s",
                               n.spde = spde$n.spde,
                               n.group = timesn)


A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(ww_all[,c('X1', 'Y1')]), group = group)

# put data into INLA stack for convinience



stk.e = 
  inla.stack(
    tag = "est",
    data = list(y = ww_all$vl_stand),
    A = list(1,A),
    effects = list(data.frame(beta = rep(1, length(ww_all$ara_id)),
                              below_loq=as.numeric(ww_all$below_loq), 
                              below_lod=as.numeric(ww_all$below_lod), 
                              lab_method=ww_all$lab_method,
                              day=ww_all$day, 
                              weekend=ww_all$weekend,
                              ara1=ww_all$ara1, 
                              day1=ww_all$day1, 
                              ara2=ww_all$ara2), 
                   s = indexs)
  )




#rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.9)))


# model formula 

formula = y ~ 0 + beta +  
  f(below_loq,model="iid") +
  f(below_lod,model="iid") +
  f(day,model="rw2", scale.model=TRUE, constr=TRUE,
    hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
  f(hol,model="linear",mean.linear=0,prec.linear=.2) +
  f(s, 
    model = spde, 
    group = s.group, 
    control.group = list(model = "ar1", hyper = list(theta = list(prior = "pccor1", param = c(0, 0.9)))))



#' ## Model A5.3.1: unique national trend
#' 
#' We directly apply model A5.3 to all ARAs. We assume one temporal trend for Switzerland, and each ARA is free to deviate from it independently.
#' 
#+ ma5.3.1, fig.width=8, fig.height=12,  R.options = list(width = 1000)
if(controls$rerun_models) {
  ms5.3.2 = inla(formula,
                 family='gaussian',
                 data = inla.stack.data(stk.e),
                 control.compute = list(config=TRUE),
                 control.predictor = list(
                   compute = TRUE,
                   A = inla.stack.A(stk.e)
                 )
                 
  )
  saveRDS(ms5.3.1,file=paste0("../",controls$savepoint,"ma5.3.1.rds"))
} else {
  ms5.3.1 = readRDS(file=paste0("../",controls$savepoint,"ma5.3.1.rds"))
}

dp_s = copy(ww_all)

index_s <- inla.stack.index(stack = stk.e, tag = "est")$data
dp_s$pred_mean = ms5.3.2$summary.fitted.values[index_s, "mean"]
dp_s$pred_lower = ms5.3.2$summary.fitted.values[index_s, "0.025quant"]
dp_s$pred_upper = ms5.3.2$summary.fitted.values[index_s, "0.975quant"]


index_e <- inla.stack.index(stack = stk.e, tag = "est")$effects


ms5.3.2$summary.fitted.values[index_e,]


ww_all = cbind(ww_all, ms5.3.2$summary.fitted.values[index_e,])

ggplot(dp_s[below_lod==0 & below_loq==0, ]) + 
  geom_point(aes(x=day1, y=vl_stand), size=0.3) + 
  geom_line(aes(x=day1, y=pred_mean), color='red', linewidth=0.2, alpha=0.8)+
  geom_ribbon(aes(x=day1, ymin=pred_lower, ymax=pred_upper), fill='red', alpha=0.2)+
  facet_wrap(~ara_id)+
  scale_y_continuous(trans='log')

# model formula 

formula = y ~ 0 + beta +  
  f(below_loq,model="iid") +
  f(below_lod,model="iid") +
  f(day,model="rw2", scale.model=TRUE, constr=TRUE,
    hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
  f(s, 
    model = spde, 
    group = s.group, 
    control.group = list(model = "rw1", hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))))



#' ## Model A5.3.1: unique national trend
#' 
#' We directly apply model A5.3 to all ARAs. We assume one temporal trend for Switzerland, and each ARA is free to deviate from it independently.
#' 
#+ ma5.3.1, fig.width=8, fig.height=12,  R.options = list(width = 1000)
if(controls$rerun_models) {
  ms5.3.2 = inla(formula,
                 family='gaussian',
                 data = inla.stack.data(stk.e),
                 control.compute = list(config=TRUE),
                 control.predictor = list(
                   compute = TRUE,
                   A = inla.stack.A(stk.e)
                 )
                 
  )
  saveRDS(ms5.3.1,file=paste0("../",controls$savepoint,"ma5.3.1.rds"))
} else {
  ms5.3.1 = readRDS(file=paste0("../",controls$savepoint,"ma5.3.1.rds"))
}

formula = y ~ 0 + beta +  
  f(below_loq,model="iid") +
  f(below_lod,model="iid") +
  f(day,model="rw2", scale.model=TRUE, constr=TRUE,
    hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
  f(s, 
    model = spde, 
    group = s.group, 
    control.group = list(model = "rw1", hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))))


if(controls$rerun_models) {
  ms5.3.4 = inla(formula,
                 family='gamma',
                 data = inla.stack.data(stk.e),
                 control.compute = list(config=TRUE),
                 control.inla=list(cmin=0),
                 control.predictor = list(
                   link=1,
                   compute = TRUE,
                   A = inla.stack.A(stk.e)
                 ), 
                 safe=TRUE
                 
  )
  saveRDS(ms5.3.1,file=paste0("../",controls$savepoint,"ma5.3.1.rds"))
} else {
  ms5.3.1 = readRDS(file=paste0("../",controls$savepoint,"ma5.3.1.rds"))
}


dp_s = copy(ww_all)

index_s <- inla.stack.index(stack = stk.e, tag = "est")$data
dp_s$pred_mean = ms5.3.4$summary.fitted.values[index_s, "mean"]
dp_s$pred_lower = ms5.3.4$summary.fitted.values[index_s, "0.025quant"]
dp_s$pred_upper = ms5.3.4$summary.fitted.values[index_s, "0.975quant"]


ms4_plot = ggplot(dp_s[below_lod==0 & below_loq==0, ]) + 
  geom_point(aes(x=day1, y=vl_stand), size=0.3) + 
  geom_line(aes(x=day1, y=pred_mean), color='red', linewidth=0.2, alpha=0.8)+
  geom_ribbon(aes(x=day1, ymin=pred_lower, ymax=pred_upper), fill='red', alpha=0.2)+
  facet_wrap(~ara_id)+
  scale_y_continuous(trans='log')
