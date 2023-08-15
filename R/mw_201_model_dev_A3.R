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

# correspondence table
corr_all = ww_all %>% 
  group_by(ara_n,ara_id,ara_name,kt,pop,lab,lab_n,lab2,lab_n2,lab_method,lab_method_n,ara1,ara2,NUTS2_name) %>% 
  count() %>% 
  ungroup() 
corr_all_ara = ww_all %>% 
  group_by(ara1,ara_name,ara_id,kt,NUTS2_name) %>% 
  count() %>% 
  ungroup() 
saveRDS(corr_all_ara,file=paste0("../",controls$savepoint,"corr_all_ara.rds"))

if(!controls$rerun_models) {
  mw_100_desc_table(ww_all) %>%
    dplyr::mutate(across(everything(),as.character)) %>%
    tidyr::gather() %>%
    flextable::flextable(cwidth=c(4,4))
}

#' ## Model A5.3.1: unique national trend
#' 
#' We directly apply model A5.3 to all ARAs. We assume one temporal trend for Switzerland, and each ARA is free to deviate from it independently.
#' 
#+ ma5.3.1, fig.width=8, fig.height=12,  R.options = list(width = 1000)
if(controls$rerun_models) {
  ma5.3.1 = INLA::inla(vl ~ 1 +
                         f(below_loq,model="iid") +
                         f(below_lod,model="iid") +
                         f(day,model="rw2", scale.model=TRUE, constr=TRUE,
                           hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                         f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
                         f(hol,model="linear",mean.linear=0,prec.linear=.2) +
                         f(ara1,model="iid") +
                         f(day1,model="rw1", scale.model=TRUE, constr=TRUE,
                           group=ara2, control.group=list(model="iid"),
                           hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))),
                       data = ww_all,
                       family = "gamma",
                       control.compute = list(waic=TRUE,config=TRUE),
                       control.predictor = list(compute=TRUE,link=1))
  saveRDS(ma5.3.1,file=paste0("../",controls$savepoint,"ma5.3.1.rds"))
} else {
  ma5.3.1 = readRDS(file=paste0("../",controls$savepoint,"ma5.3.1.rds"))
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

#' 
#' ## Model A5.3.2: effect of lab and method change
#' 
#' We add a covariate to measure the effect of the lab and of changes in methodology. To allow identifiability, we have to make sure that there are multiple ARAs per laboratory, so we group together KLBS (only 1 ARA) and KLZH (12 ARAs). The reference lab is ALTGR.
#'  
#+ ma5.3.2a, fig.width=8, fig.height=12,  R.options = list(width = 1000)
if(controls$rerun_models) {
  ma5.3.2 = INLA::inla(vl ~ 1 +
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
                         lab_method,
                       data = ww_all,
                       family = "gamma",
                       control.compute = list(waic=TRUE,config=TRUE),
                       control.predictor = list(compute=TRUE,link=1))
  saveRDS(ma5.3.2,file=paste0("../",controls$savepoint,"ma5.3.2.rds"))
} else {
  ma5.3.2 = readRDS(file=paste0("../",controls$savepoint,"ma5.3.2.rds"))
}
summary(ma5.3.2)
summary_exp_vl(ma5.3.2,pars="lab|method|hol|weekend")
ppp_vl_ara(ww_all,ma5.3.2)
#+ ma5.3.2b, fig.width=8, fig.height=4,  R.options = list(width = 1000)
avg_time_trend(ww_all,ma5.3.2)
#+ ma5.3.2c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_130_map_relative_vl(ma5.3.2,corr_all_ara,shapes)
#+ ma5.3.2d, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_131_map_deviation_from_average(ma5.3.2,corr_all_ara,ww_all,shapes,12)

#' 
#' ## Model A5.3.3: regional effects
#' 
#' Instead of just one temporal trend for Switzerland, we now allow independent temporal trends for each NUTS-2 region. ARAs- within each region are then allowed to deviate from the regional trend independently.
#' 
#+ ma5.3.3a, fig.width=8, fig.height=12,  R.options = list(width = 1000)
if(controls$rerun_models) {
  ma5.3.3 = INLA::inla(vl ~ 1 +
                         f(below_loq,model="iid") +
                         f(below_lod,model="iid") +
                         f(day,model="rw2", scale.model=TRUE, constr=TRUE,
                           group=NUTS2, control.group=list(model="iid"),
                           hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                         f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
                         f(hol,model="linear",mean.linear=0,prec.linear=.2) +
                         f(ara1,model="iid") +
                         f(day1,model="rw1", scale.model=TRUE, constr=TRUE,
                           group=ara2, control.group=list(model="iid"),
                           hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                         lab_method,
                       data = ww_all,
                       family = "gamma",
                       control.compute = list(waic=TRUE,config=TRUE),
                       control.predictor = list(compute=TRUE,link=1))
  saveRDS(ma5.3.3,file=paste0("../",controls$savepoint,"ma5.3.3.rds"))
} else {
  ma5.3.3 = readRDS(file=paste0("../",controls$savepoint,"ma5.3.3.rds"))
}
summary(ma5.3.3)
summary_exp_vl(ma5.3.3,pars="lab|method|hol|weekend")
ppp_vl_ara(ww_all,ma5.3.3)
#+ ma5.3.3b, fig.width=8, fig.height=4,  R.options = list(width = 1000)
avg_time_trend_reg(ww_all,ma5.3.3)
#+ ma5.3.3c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_130_map_relative_vl(ma5.3.3,corr_all_ara,shapes)
#+ ma5.3.3d, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_131_map_deviation_from_average(ma5.3.3,corr_all_ara,ww_all,shapes,12)


#' 
#' ## Model A5.3.3b: new ARAs
#' 
#' 
#' 
#+ ma5.3.3a, fig.width=8, fig.height=12,  R.options = list(width = 1000)
if(controls$rerun_models) {
  ma5.3.3 = INLA::inla(vl ~ 1 +
                         f(below_loq,model="iid") +
                         f(below_lod,model="iid") +
                         f(day,model="rw2", scale.model=TRUE, constr=TRUE,
                           group=NUTS2, control.group=list(model="iid"),
                           hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                         f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
                         f(hol,model="linear",mean.linear=0,prec.linear=.2) +
                         f(ara1,model="iid") +
                         f(day1,model="rw1", scale.model=TRUE, constr=TRUE,
                           group=ara2, control.group=list(model="iid"),
                           hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                         lab_method,
                       data = ww_all,
                       family = "gamma",
                       control.compute = list(waic=TRUE,config=TRUE),
                       control.predictor = list(compute=TRUE,link=1))
  saveRDS(ma5.3.3,file=paste0("../",controls$savepoint,"ma5.3.3.rds"))
} else {
  ma5.3.3 = readRDS(file=paste0("../",controls$savepoint,"ma5.3.3.rds"))
}
summary(ma5.3.3)
summary_exp_vl(ma5.3.3,pars="lab|method|hol|weekend")
ppp_vl_ara(ww_all,ma5.3.3)
#+ ma5.3.3b, fig.width=8, fig.height=4,  R.options = list(width = 1000)
avg_time_trend_reg(ww_all,ma5.3.3)
#+ ma5.3.3c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_130_map_relative_vl(ma5.3.3,corr_all_ara,shapes)
#+ ma5.3.3d, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_131_map_deviation_from_average(ma5.3.3,corr_all_ara,ww_all,shapes,12)


#' 
#' 
#' #' ## Model A5.4.1: spatial structure
#' #' 
#' #' 
#' #+ ma5.4.1a, fig.width=8, fig.height=12,  R.options = list(width = 1000)
#' if(controls$rerun_models) {
#'   shapes_all = shapes$ara_shp %>% 
#'     dplyr::filter(ara_id %in% ww_all$ara_id) %>% 
#'     left_join(corr_all_ara,by = join_by(ara_id)) %>% 
#'     dplyr::arrange(ara1)
#'   sf_use_s2(FALSE)
#'   graph_all = spdep::poly2nb(shapes_all)
#'   path_graph = paste0("../",controls$savepoint,"W_all.adj")
#'   nb2INLA(path_graph, graph_all)
#'   ma5.4.1 = INLA::inla(vl ~ 1 +
#'                          f(below_loq,model="iid") +
#'                          f(below_lod,model="iid") +
#'                          f(day,model="rw2", scale.model=TRUE, constr=TRUE,
#'                            group=NUTS2, control.group=list(model="iid"),
#'                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#'                          f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
#'                          f(hol,model="linear",mean.linear=0,prec.linear=.2) +
#'                          f(ara1,model="bym2",
#'                            graph=path_graph,
#'                            scale.model = TRUE, constr = TRUE, 
#'                            hyper = list(theta1 = list("PCprior", c(1, 0.01)), 
#'                                         theta2 = list("PCprior", c(0.5, 0.5)))) +
#'                          f(day1,model="rw1", scale.model=TRUE, constr=TRUE,
#'                            group=ara2, control.group=list(model="iid"),
#'                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#'                          lab_method,
#'                        data = ww_all,
#'                        family = "gamma",
#'                        control.compute = list(waic=TRUE,config=TRUE),
#'                        control.predictor = list(compute=TRUE,link=1))
#'   saveRDS(ma5.4.1,file=paste0("../",controls$savepoint,"ma5.4.1.rds"))
#' } else {
#'   ma5.4.1 = readRDS(file=paste0("../",controls$savepoint,"ma5.4.1.rds"))
#' }
#' summary(ma5.4.1)
#' summary_exp_vl(ma5.4.1,pars="lab|method|hol|weekend")
#' ppp_vl_ara(ww_all,ma5.4.1)
#' #+ ma5.4.1b, fig.width=8, fig.height=6,  R.options = list(width = 1000)
#' mw_130_map_relative_vl(ma5.4.1,corr_all_ara,shapes)
#' #+ ma5.4.1c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
#' mw_131_map_deviation_from_average(ma5.4.1,corr_all_ara,ww_all,shapes,10)
