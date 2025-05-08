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
if(!exists("controls")) controls = readRDS(fs::path("../savepoints/savepoint_2025-01-24/controls.rds"))
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
  dplyr::mutate(lab_method=factor(paste0(lab2,"_",method)),
                lab_method=relevel(lab_method, ref="EAWAG_0"),
                lab_method_n=as.numeric(lab_method),
                prop_under_20b=(prop_under_20-mean(prop_under_20))/sd(prop_under_20),
                prop_over_65b=(prop_over_65-mean(prop_over_65))/sd(prop_over_65),
                prop_non_ch_eub=(prop_non_ch_eu-mean(prop_non_ch_eu))/sd(prop_non_ch_eu),
                employment_factorb=(employment_factor-mean(employment_factor))/sd(employment_factor),
                ssep3_medb=(ssep3_med-mean(ssep3_med))/sd(ssep3_med))
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


#' ## Model A5.4.0: baseline
#' 
#' We directly apply model A5.4 to all ARAs. We assume one temporal trend for Switzerland, and each ARA is free to deviate from it independently.
#' 
#+ ma5.4.0, fig.width=8, fig.height=12,  R.options = list(width = 1000)
# if(FALSE) {
#   ma5.4.0 = INLA::inla(vl ~ 1 +
#                          f(below_loq,model="iid") +
#                          f(below_lod,model="iid") +
#                          f(day,model="rw2", scale.model=TRUE, constr=TRUE,
#                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#                          f(ara1,model="iid") +
#                          f(day1,model="rw1", scale.model=TRUE, constr=TRUE,
#                            group=ara2, control.group=list(model="iid"),
#                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))),
#                        data = ww_all,
#                        family = "gamma",
#                        control.compute = list(waic=TRUE,config=TRUE),
#                        control.predictor = list(compute=TRUE,link=1),
#                        num.threads = "4:1")
#   saveRDS(ma5.4.0,file=paste0("../",controls$savepoint,"ma5.4.0.rds"))
# } else {
ma5.4.0 = readRDS(file=paste0("../",controls$savepoint,"ma5.4.0.rds"))

summary(ma5.4.0)
summary_exp_vl(ma5.4.0,pars="lab|method|hol|weekend|pop_total|prop_under_20|prop_over_65|prop_non_ch_eu|ssep3_|emp")
#+ ma5.4.0a, fig.width=8, fig.height=12,  R.options = list(width = 1000)
ppp_vl_ara(ww_all,ma5.4.0)
#+ ma5.4.0b, fig.width=8, fig.height=4,  R.options = list(width = 1000)
avg_time_trend(ww_all,ma5.4.0)
#+ ma5.4.0c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_130_map_relative_vl(ma5.4.0,corr_all_ara,shapes)
#+ ma5.4.0d, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_131_map_deviation_from_average(ma5.4.0,corr_all_ara,ww_all,shapes,12)
# }

#' ## Model A5.4.1: unique national trend
#' 
#' We directly apply model A5.4 to all ARAs. We assume one temporal trend for Switzerland, and each ARA is free to deviate from it independently.
#' 
# if(FALSE) {
#   ma5.4.1 = INLA::inla(vl ~ 1 +
#                          f(below_loq,model="iid") +
#                          f(below_lod,model="iid") +
#                          f(day,model="rw2", scale.model=TRUE, constr=TRUE,
#                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#                          f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(hol,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(ara1,model="iid") +
#                          f(day1,model="rw1", scale.model=TRUE, constr=TRUE,
#                            group=ara2, control.group=list(model="iid"),
#                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#                          f(prop_under_20b,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(prop_over_65b,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(prop_non_ch_eub,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(ssep3_medb,model="linear",mean.linear=0,prec.linear=.2)  +
#                          f(employment_factorb,model="linear",mean.linear=0,prec.linear=.2),
#                        data = ww_all,
#                        family = "gamma",
#                        control.compute = list(waic=TRUE,config=TRUE),
#                        control.predictor = list(compute=TRUE,link=1),
#                        num.threads = "4:1")
#   saveRDS(ma5.4.1,file=paste0("../",controls$savepoint,"ma5.4.1.rds"))
# } else {
ma5.4.1 = readRDS(file=paste0("../",controls$savepoint,"ma5.4.1.rds"))

summary(ma5.4.1)
summary_exp_vl(ma5.4.1,pars="lab|method|hol|weekend|pop_total|prop_under_20|prop_over_65|prop_non_ch_eu|ssep3_|emp")
#+ ma5.4.1a, fig.width=8, fig.height=12,  R.options = list(width = 1000)
ppp_vl_ara(ww_all,ma5.4.1)
#+ ma5.4.1b, fig.width=8, fig.height=4,  R.options = list(width = 1000)
avg_time_trend(ww_all,ma5.4.1)
#+ ma5.4.1c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_130_map_relative_vl(ma5.4.1,corr_all_ara,shapes)
#+ ma5.4.1d, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_131_map_deviation_from_average(ma5.4.1,corr_all_ara,ww_all,shapes,12)
# }
#' 
#' ## Model A5.4.2: effect of lab and method change
#' 
#' We add a covariate to measure the effect of the lab and of changes in methodology. 
#' To allow identifiability, we have to make sure that there are multiple ARAs per laboratory, 
#' so we group together KLBS (only 1 ARA) and KLZH (12 ARAs). The reference lab is EAWAG.
#'  
# if(FALSE) {
#   ma5.4.2 = INLA::inla(vl ~ 1 +
#                          f(below_loq,model="iid") +
#                          f(below_lod,model="iid") +
#                          f(day,model="rw2", scale.model=TRUE, constr=TRUE,
#                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#                          f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(hol,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(ara1,model="iid") +
#                          f(day1,model="rw1", scale.model=TRUE, constr=TRUE,
#                            group=ara2, control.group=list(model="iid"),
#                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#                          lab_method +
#                          f(prop_under_20b,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(prop_over_65b,model="linear",mean.linear=0,prec.linear=.2) +
#                          # f(prop_non_ch_eub,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(ssep3_medb,model="linear",mean.linear=0,prec.linear=.2)  +
#                          f(employment_factorb,model="linear",mean.linear=0,prec.linear=.2),
#                        data = ww_all,
#                        family = "gamma",
#                        control.compute = list(waic=TRUE,config=TRUE),
#                        control.predictor = list(compute=TRUE,link=1),
#                        num.threads = "4:1")
#   saveRDS(ma5.4.2,file=paste0("../",controls$savepoint,"ma5.4.2b.rds"))
# } else {
ma5.4.2 = readRDS(file=paste0("../",controls$savepoint,"ma5.4.2b.rds"))

summary(ma5.4.2)
summary_exp_vl(ma5.4.2,pars="lab|method|hol|weekend|pop_total|prop_under_20|prop_over_65|prop_non_ch_eu|ssep3_|emp")
#+ ma5.4.2a, fig.width=8, fig.height=12,  R.options = list(width = 1000)
ppp_vl_ara(ww_all,ma5.4.2)
#+ ma5.4.2b, fig.width=8, fig.height=4,  R.options = list(width = 1000)
avg_time_trend(ww_all,ma5.4.2)
#+ ma5.4.2c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_130_map_relative_vl(ma5.4.2,corr_all_ara,shapes)
#+ ma5.4.2d, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_131_map_deviation_from_average(ma5.4.2,corr_all_ara,ww_all,shapes,12)
# }
#' 
#' ## Model A5.4.3: regional effects
#' 
#' Instead of just one temporal trend for Switzerland, we now allow independent temporal trends for each NUTS-2 region. ARAs- within each region are then allowed to deviate from the regional trend independently.
#' 
# if(FALSE) {
#   ma5.4.3 = INLA::inla(vl ~ 1 +
#                          f(below_loq,model="iid") +
#                          f(below_lod,model="iid") +
#                          f(day,model="rw2", scale.model=TRUE, constr=TRUE,
#                            group=NUTS2, control.group=list(model="iid"),
#                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#                          f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(hol,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(ara1,model="iid") +
#                          f(day1,model="rw1", scale.model=TRUE, constr=TRUE,
#                            group=ara2, control.group=list(model="iid"),
#                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01))))  +
#                          lab_method +
#                          f(prop_under_20b,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(prop_over_65b,model="linear",mean.linear=0,prec.linear=.2) +
#                          # f(prop_non_ch_eub,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(ssep3_medb,model="linear",mean.linear=0,prec.linear=.2)  +
#                          f(employment_factorb,model="linear",mean.linear=0,prec.linear=.2),
#                        data = ww_all,
#                        family = "gamma",
#                        control.compute = list(waic=TRUE,config=TRUE),
#                        control.predictor = list(compute=TRUE,link=1),
#                        num.threads="4:1")
#   saveRDS(ma5.4.3,file=paste0("../",controls$savepoint,"ma5.4.3.rds"))
# } else {
ma5.4.3 = readRDS(file=paste0("../",controls$savepoint,"ma5.4.3.rds"))

summary(ma5.4.3)
summary_exp_vl(ma5.4.3,pars="lab|method|hol|weekend|pop_total|prop_under_20|prop_over_65|prop_non_ch_eu|ssep3_|emp")
#+ ma5.4.3a, fig.width=8, fig.height=12,  R.options = list(width = 1000)
ppp_vl_ara(ww_all,ma5.4.3)
#+ ma5.4.3b, fig.width=8, fig.height=4,  R.options = list(width = 1000)
avg_time_trend_reg(ww_all,ma5.4.3)
#+ ma5.4.3c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_130_map_relative_vl(ma5.4.3,corr_all_ara,shapes)
#+ ma5.4.3d, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_131_map_deviation_from_average(ma5.4.3,corr_all_ara,ww_all,shapes,12)
# }

#' 
#' ## Model A5.4.4: spatial structure
#' 
#' Instead of just one temporal trend for Switzerland, we now allow independent temporal trends for each NUTS-2 region. ARAs- within each region are then allowed to deviate from the regional trend independently.
#' 
# if(FALSE) {
#   shapes_all = shapes$ara_shp %>%
#     dplyr::filter(ara_id %in% ww_all$ara_id) %>%
#     left_join(corr_all_ara,by = join_by(ara_id)) %>%
#     dplyr::arrange(ara1)
#   sf_use_s2(FALSE)
#   graph_all = spdep::poly2nb(shapes_all)
#   path_graph = paste0("../",controls$savepoint,"W_all.adj")
#   nb2INLA(path_graph, graph_all)
#   if(FALSE) {
#     plot(st_geometry(shapes_all), border="grey")
#     coords = st_coordinates(st_centroid(st_geometry(shapes_all)))
#     plot(graph_all, coords, add=TRUE, col="red")
#   }
#   
#   ma5.4.4 = INLA::inla(vl ~ 1 +
#                          f(below_loq,model="iid") +
#                          f(below_lod,model="iid") +
#                          f(day,model="rw2", scale.model=TRUE, constr=TRUE,
#                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#                          f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(hol,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(ara1,model="bym2",
#                            graph=path_graph,
#                            scale.model = TRUE, constr = TRUE,
#                            hyper = list(theta1 = list("PCprior", c(1, 0.01)),
#                                         theta2 = list("PCprior", c(0.5, 0.5)))) +
#                          f(day1,model="rw1", scale.model=TRUE, constr=TRUE,
#                            group=ara2, control.group=list(model="iid"),
#                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01))))  +
#                          lab_method +
#                          f(prop_under_20b,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(prop_over_65b,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(ssep3_medb,model="linear",mean.linear=0,prec.linear=.2)  +
#                          f(employment_factorb,model="linear",mean.linear=0,prec.linear=.2),
#                        data = ww_all,
#                        family = "gamma",
#                        control.inla=list(cmin=0),
#                        control.compute = list(waic=TRUE,config=TRUE),
#                        control.predictor = list(compute=TRUE,link=1),
#                        num.threads="4:1")
#   saveRDS(ma5.4.4,file=paste0("../",controls$savepoint,"ma5.4.4.rds"))
# } else {
ma5.4.4 = readRDS(file=paste0("../",controls$savepoint,"ma5.4.4.rds"))

summary(ma5.4.4)
summary_exp_vl(ma5.4.4,pars="lab|method|hol|weekend|pop_total|prop_under_20|prop_over_65|prop_non_ch_eu|ssep3_|emp")
#+ ma5.4.4a, fig.width=8, fig.height=12,  R.options = list(width = 1000)
ppp_vl_ara(ww_all,ma5.4.4)
#+ ma5.4.4b, fig.width=8, fig.height=4,  R.options = list(width = 1000)
avg_time_trend(ww_all,ma5.4.4)
#+ ma5.4.4c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_130_map_relative_vl(ma5.4.4,corr_all_ara,shapes)
#+ ma5.4.4d, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_131_map_deviation_from_average(ma5.4.4,corr_all_ara,ww_all,shapes,12)
# }

#' 
#' ## Model A5.4.5: spatial structure + regional effects
#' 
#' Instead of just one temporal trend for Switzerland, we now allow independent temporal trends for each NUTS-2 region. ARAs- within each region are then allowed to deviate from the regional trend independently.
#' 
# if(FALSE) {
#   ma5.4.5 = INLA::inla(vl ~ 1 +
#                          f(below_loq,model="iid") +
#                          f(below_lod,model="iid") +
#                          f(day,model="rw2", scale.model=TRUE, constr=TRUE,
#                            group=NUTS2, control.group=list(model="iid"),
#                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#                          f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(hol,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(ara1,model="bym2",
#                            graph=path_graph,
#                            scale.model = TRUE, constr = TRUE,
#                            hyper = list(theta1 = list("PCprior", c(1, 0.01)),
#                                         theta2 = list("PCprior", c(0.5, 0.5)))) +
#                          f(day1,model="rw1", scale.model=TRUE, constr=TRUE,
#                            group=ara2, control.group=list(model="iid"),
#                            hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01))))  +
#                          lab_method +
#                          f(prop_under_20b,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(prop_over_65b,model="linear",mean.linear=0,prec.linear=.2) +
#                          f(ssep3_medb,model="linear",mean.linear=0,prec.linear=.2)  +
#                          f(employment_factorb,model="linear",mean.linear=0,prec.linear=.2),
#                        data = ww_all,
#                        family = "gamma",
#                        control.inla=list(cmin=0),
#                        control.compute = list(waic=TRUE,config=TRUE),
#                        control.predictor = list(compute=TRUE,link=1),
#                        num.threads="4:1")
#   saveRDS(ma5.4.5,file=paste0("../",controls$savepoint,"ma5.4.5.rds"))
# } else {
ma5.4.5 = readRDS(file=paste0("../",controls$savepoint,"ma5.4.5.rds"))

summary(ma5.4.5)
summary_exp_vl(ma5.4.5,pars="lab|method|hol|weekend|pop_total|prop_under_20|prop_over_65|prop_non_ch_eu|ssep3_|emp")
#+ ma5.4.5a, fig.width=8, fig.height=12
ppp_vl_ara(ww_all,ma5.4.5)
#+ ma5.4.5b, fig.width=8, fig.height=4
avg_time_trend_reg(ww_all,ma5.4.5)
#+ ma5.4.5c, fig.width=8, fig.height=6
mw_130_map_relative_vl(ma5.4.5,corr_all_ara,shapes)
#+ ma5.4.5d, fig.width=8, fig.height=6
mw_131_map_deviation_from_average(ma5.4.5,corr_all_ara,ww_all,shapes,12)
# }

