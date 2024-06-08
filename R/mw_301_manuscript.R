#' ---
#' title: "Tables and figures for manuscript"
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
source("setup.R")
ww1 = readRDS(fs::path("../",controls$savepoint,"ww1.rds"))
shapes = readRDS(fs::path("../",controls$savepoint,"shapes.rds"))
ww_all = readRDS(file=paste0("../",controls$savepoint,"ww_all.rds"))

#' 
#' # Description
#' 

mw_100_desc_table(ww_all) %>%
  dplyr::mutate(across(everything(),as.character)) %>%
  tidyr::gather() %>%
  xtable::xtable() %>% 
  print(include.rownames = FALSE)


ww1 %>% 
  group_by(ara_name) %>% 
  summarise(n=n(),
            min=min(date),
            max=max(date),
            dur=as.numeric(max-min),
            freq=n/dur*7) %>% view()

ww1$kt %>% unique() %>% length()

#' # Figure 1
#' 
#' Version 1: separate counts and VL
g1 = mw_110_map_missing(ww1,shapes)
g2 = mw_106_fig_vl_time(ww1)

cowplot::plot_grid(g1,g2,labels=c("A","B"))
ggsave(file.path("..",controls$savepoint,"figures_paper","fig1.pdf"),width=8,height=3)

#' Version 2: together but median by week
mw_104_fig_vl(ww1)
ggsave(file.path("..",controls$savepoint,"figures_paper","fig1b.pdf"),width=6,height=8)

#' ## Model A5.4.1: unique national trend
#' 
#' We directly apply model A5.4 to all ARAs. We assume one temporal trend for Switzerland, and each ARA is free to deviate from it independently.
#' 
#+ ma5.4.1, fig.width=8, fig.height=12,  R.options = list(width = 1000)
if(controls$rerun_models) {
  ma5.4.1 = INLA::inla(vl ~ 1 +
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
                         f(prop_under_20b,model="linear",mean.linear=0,prec.linear=.2) +
                         f(prop_over_65b,model="linear",mean.linear=0,prec.linear=.2) +
                         f(prop_non_ch_eub,model="linear",mean.linear=0,prec.linear=.2) +
                         f(ssep3_med,model="linear",mean.linear=0,prec.linear=.2),
                       data = ww_all,
                       family = "gamma",
                       control.compute = list(waic=TRUE,config=TRUE),
                       control.predictor = list(compute=TRUE,link=1),
                       num.threads = "4:1")
  saveRDS(ma5.4.1,file=paste0("../",controls$savepoint,"ma5.4.1.rds"))
} else {
  ma5.4.1 = readRDS(file=paste0("../",controls$savepoint,"ma5.4.1.rds"))
}
summary(ma5.4.1)
summary_exp_vl(ma5.4.1,pars="lab|method|hol|weekend|pop_total|prop_under_20|prop_over_65|prop_non_ch_eu|ssep3_")
ppp_vl_ara(ww_all,ma5.4.1)
#+ ma5.4.1b, fig.width=8, fig.height=4,  R.options = list(width = 1000)
avg_time_trend(ww_all,ma5.4.1)
#+ ma5.4.1c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_130_map_relative_vl(ma5.4.1,corr_all_ara,shapes)
#+ ma5.4.1d, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_131_map_deviation_from_average(ma5.4.1,corr_all_ara,ww_all,shapes,12)

#' 
#' ## Model A5.4.2: effect of lab and method change
#' 
#' We add a covariate to measure the effect of the lab and of changes in methodology. 
#' To allow identifiability, we have to make sure that there are multiple ARAs per laboratory, 
#' so we group together KLBS (only 1 ARA) and KLZH (12 ARAs). The reference lab is EAWAG.
#'  
#+ ma5.4.2a, fig.width=8, fig.height=12,  R.options = list(width = 1000)
if(controls$rerun_models) {
  ma5.4.2 = INLA::inla(vl ~ 1 +
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
                         f(lab_method,model="linear",mean.linear=0,prec.linear=.2) +
                         f(prop_under_20b,model="linear",mean.linear=0,prec.linear=.2) +
                         f(prop_over_65b,model="linear",mean.linear=0,prec.linear=.2) +
                         f(prop_non_ch_eub,model="linear",mean.linear=0,prec.linear=.2) +
                         f(ssep3_med,model="linear",mean.linear=0,prec.linear=.2),
                       data = ww_all,
                       family = "gamma",
                       control.compute = list(waic=TRUE,config=TRUE),
                       control.predictor = list(compute=TRUE,link=1),
                       num.threads = "4:1")
  saveRDS(ma5.4.2,file=paste0("../",controls$savepoint,"ma5.4.2.rds"))
} else {
  ma5.4.2 = readRDS(file=paste0("../",controls$savepoint,"ma5.4.2.rds"))
}
summary(ma5.4.2)
summary_exp_vl(ma5.4.2,pars="lab|method|hol|weekend|pop_total|prop_under_20|prop_over_65|prop_non_ch_eu|ssep3_")
ppp_vl_ara(ww_all,ma5.4.2)
#+ ma5.4.2b, fig.width=8, fig.height=4,  R.options = list(width = 1000)
avg_time_trend(ww_all,ma5.4.2)
#+ ma5.4.2c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_130_map_relative_vl(ma5.4.2,corr_all_ara,shapes)
#+ ma5.4.2d, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_131_map_deviation_from_average(ma5.4.2,corr_all_ara,ww_all,shapes,12)

#' 
#' ## Model A5.4.3: regional effects
#' 
#' Instead of just one temporal trend for Switzerland, we now allow independent temporal trends for each NUTS-2 region. ARAs- within each region are then allowed to deviate from the regional trend independently.
#' 
#+ ma5.4.3a, fig.width=8, fig.height=12,  R.options = list(width = 1000)
if(controls$rerun_models) {
  ma5.4.3 = INLA::inla(vl ~ 1 +
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
                           hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01))))  +
                         f(lab_method,model="linear",mean.linear=0,prec.linear=.2) +
                         f(prop_under_20b,model="linear",mean.linear=0,prec.linear=.2) +
                         f(prop_over_65b,model="linear",mean.linear=0,prec.linear=.2) +
                         f(prop_non_ch_eub,model="linear",mean.linear=0,prec.linear=.2) +
                         f(ssep3_med,model="linear",mean.linear=0,prec.linear=.2),
                       data = ww_all,
                       family = "gamma",
                       control.compute = list(waic=TRUE,config=TRUE),
                       control.predictor = list(compute=TRUE,link=1),
                       num.threads="4:1")
  saveRDS(ma5.4.3,file=paste0("../",controls$savepoint,"ma5.4.3.rds"))
} else {
  ma5.4.3 = readRDS(file=paste0("../",controls$savepoint,"ma5.4.3.rds"))
}
summary(ma5.4.3)
summary_exp_vl(ma5.4.3,pars="lab|method|hol|weekend|pop_total|prop_under_20|prop_over_65|prop_non_ch_eu|ssep3_")
ppp_vl_ara(ww_all,ma5.4.3)
#+ ma5.4.3b, fig.width=8, fig.height=4,  R.options = list(width = 1000)
avg_time_trend_reg(ww_all,ma5.4.3)
#+ ma5.4.3c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_130_map_relative_vl(ma5.4.3,corr_all_ara,shapes)
#+ ma5.4.3d, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_131_map_deviation_from_average(ma5.4.3,corr_all_ara,ww_all,shapes,12)



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
