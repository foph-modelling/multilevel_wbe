#' ---
#' title: "WBE for SARS-CoV-2 in Switzerland: model development (A4) "
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

#' 
#' The objective of this section is to explore using log-Gaussian Cox-processes (LGCP) instead of BYM to model the geographical structure.
#' The point is to use the actual distance between the centroïds of each ARAs instead of the neighboring structure. 
#' For this, we follow the same structure as the A2 section, and focus on one region.
#' 
#' # Mittelland
#' 
#' 

# select NUTS2 region
select_NUTS2 = "Mittelland"
ww_reg = ww1 %>%
  # log
  dplyr::mutate(logvl=log(vl)) %>%
  mutate(vl=if_else(vl==0, 1, vl)) %>%
  # select one ARA per NUTS-2
  dplyr::filter(NUTS2_name==select_NUTS2) %>%
  # create indexes for INLA
  dplyr::mutate(day1=day,
                ara1=as.numeric(as.factor(ara_n)),
                ara2=ara1)

# correspondence table
corr_reg = ww_reg %>% 
  group_by(ara_n,ara_id,ara_name,kt,pop,lab,lab_n,ara1,ara2) %>% 
  count() %>% 
  ungroup()

mw_100_desc_table(ww_reg) %>% 
  dplyr::mutate(across(everything(),as.character)) %>% 
  tidyr::gather() %>% 
  flextable::flextable(cwidth=c(4,4)) 

mw_100_desc_table(ww_reg,ara_name) %>% 
  dplyr::select(1:7)  %>% 
  flextable::flextable(cwidth=rep(4,7)) 

#+ desc_missing, fig.width=6, fig.height=3.5
mw_110_map_missing(ww_reg,shapes) + labs(title="Number of measurements")

#+ desc_vl, fig.width=8, fig.height=3.5
ggplot(ww_reg) + geom_line(aes(x=date,y=logvl,colour=ara_name),alpha=.6) + scale_colour_discrete(guide="none") + labs(title="Viral load measurements")


#' ## Model A5.3: space-time interaction
#' 
#' We restart from the A5.3 model (see section A2), with a random intercept and different time trends by ARA.
#+ ma5.3, fig.width=8, fig.height=8,  R.options = list(width = 1000)
if(controls$rerun_models) {
  ma5.3 = INLA::inla(vl ~ 1 +
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
                     data = ww_reg,
                     family = "gamma",
                     control.compute = list(waic=TRUE,config=TRUE),
                     control.predictor = list(compute=TRUE,link=1))
  saveRDS(ma5.3,file=paste0("../",controls$savepoint,"ma5.3.rds"))
} else {
  ma5.3 = readRDS(file=paste0("../",controls$savepoint,"ma5.3.rds"))
}
summary(ma5.3)
summary_exp_vl(ma5.3,pars="lab|method|hol|weekend")
ppp_vl_ara(ww_reg,ma5.3)
avg_time_trend(ww_reg,ma5.3)
print("Relative viral load by ARA compared to average:")
ma5.3$summary.random$ara1 %>% 
  as_tibble() %>% 
  dplyr::transmute(ara1=ID,
                   `exp(beta)`=round(exp(mean),2),
                   `0.025quant`=round(exp(`0.025quant`),2),
                   `0.975quant`=format(round(exp(`0.975quant`),2),scientific=FALSE)) %>% 
  dplyr::left_join(select(corr_reg,ara1,ara_name),by = join_by(ara1)) %>% 
  dplyr::arrange(-`exp(beta)`) %>% 
  select(-ara1) %>% 
  column_to_rownames(var="ara_name")
ndays = length(unique(ww_reg$day))
nara = length(unique(ww_reg$ara1))
ma5.3$summary.random$day1 %>% 
  bind_cols(day=rep(0:(ndays-1),nara),
            ara1=rep(1:nara,each=ndays)) %>% 
  left_join(corr_reg,by = join_by(ara1)) %>% 
  ggplot() +
  geom_hline(yintercept=1,linetype=2,alpha=.5) +
  geom_ribbon(aes(x=day,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=ara_name),alpha=.5) +
  geom_line(aes(x=day,y=exp(mean),colour=ara_name)) +
  facet_wrap(~ara_name) +
  scale_colour_discrete(guide="none") +
  scale_fill_discrete(guide="none") +
  scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
  coord_cartesian(ylim=c(.05,20)) +
  labs(title="Deviations from average time trend by ARA",x="Day",y="Relative viral load by ARA") 

#' 
#' ## Model A5.5: spatial structure by distance
#'  
#' We now consider the distance between the centroïds of ARAs instead of their neighbourhood structure (see https://gkonstantinoudis.github.io/INLA/T1Dmarkdown.html).
#' 
#+ fig.width=3, fig.height=3,  R.options = list(width = 1000)
# filter ARA shapes matrix
shapes_reg = shapes$ara_shp %>% 
  dplyr::filter(ara_id %in% ww_reg$ara_id) %>% 
  dplyr::left_join(corr_reg,by = join_by(ara_id)) %>% 
  dplyr::arrange(ara1)
# find centroids
centroid_reg = shapes_reg %>% 
  sf::st_centroid() %>% 
  dplyr::select(ara1,geometry)
ggplot() +
  geom_sf(data=shapes_reg,colour="black") +
  geom_sf(data=centroid_reg,colour="red")
# compute mesh
#+ fig.width=3, fig.height=3,  R.options = list(width = 1000)
max.edge = diff(range(st_coordinates(centroid_reg)[,1]))/10
mesh1 = inla.mesh.2d(loc=st_coordinates(centroid_reg),
                     max.edge =  c(1,2)*max.edge)
ggplot() +
  geom_sf(data=shapes_reg,colour="black") +
  geom_sf(data=centroid_reg,colour="red") +
  gg(mesh1) 
# create latent field
spde_reg = inla.spde2.pcmatern(mesh=mesh1, alpha=2, 
                               prior.range=c(60000,0.5), # Pr(phi < 60000) = 0.5
                               prior.sigma=c(1,0.01) # Pr(sigma > 1) = 0.01
)

#+ ma5.5, fig.width=8, fig.height=8,  R.options = list(width = 1000)
if(controls$rerun_models) {
  ma5.5 = INLA::inla(vl ~ 1 + 
                       f(day,model="rw2", scale.model=TRUE, constr=TRUE,
                        hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                       f(ara1, model = spde_reg), 
                     family='gamma', 
                     data = ww_reg,
                     control.compute=list(dic=TRUE,cpo=TRUE, config = TRUE), 
                     control.inla=list(strategy="simplified.laplace", 
                                       int.strategy="eb"))
  # 
  # ma5.5 = INLA::inla(vl ~ 1 +
  #                      f(below_loq,model="iid") +
  #                      f(below_lod,model="iid") +
  #                      f(day,model="rw2", scale.model=TRUE, constr=TRUE,
  #                        hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  #                      f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
  #                      f(hol,model="linear",mean.linear=0,prec.linear=.2) +
  #                      f(ara1,model=spde_reg) +
  #                      f(day1,model="rw1", scale.model=TRUE, constr=TRUE,
  #                        group=ara2, control.group=list(model="iid"),
  #                        hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))),
  #                    data = ww_reg,
  #                    family = "gamma",
  #                    control.compute = list(waic=TRUE,config=TRUE),
  #                    control.predictor = list(compute=TRUE,link=1))
  saveRDS(ma5.5,file=paste0("../",controls$savepoint,"ma5.5.rds"))
} else {
  ma5.5 = readRDS(file=paste0("../",controls$savepoint,"ma5.5.rds"))
}
summary(ma5.5)
summary_exp_vl(ma5.5,pars="lab|method|hol|weekend")
ppp_vl_ara(ww_reg,ma5.5) 
ma5.5$summary.random$ara1 %>% 
  as_tibble() %>% 
  dplyr::transmute(ara1=ID,
                   `exp(beta)`=round(exp(mean),2),
                   `0.025quant`=round(exp(`0.025quant`),2),
                   `0.975quant`=format(round(exp(`0.975quant`),2),scientific=FALSE)) %>% 
  dplyr::left_join(select(corr_reg,ara1,ara_name),by = join_by(ara1)) %>% 
  dplyr::arrange(-`exp(beta)`) %>% 
  dplyr::select(-ara1) %>% 
  dplyr::filter(!is.na(ara_name)) %>% 
  column_to_rownames(var="ara_name")
ndays = length(unique(ww_reg$day))
nara = length(unique(ww_reg$ara1))
ma5.5$summary.random$day1 %>% 
  bind_cols(day=rep(0:(ndays-1),nara),
            ara1=rep(1:nara,each=ndays)) %>% 
  left_join(corr_reg,by = join_by(ara1)) %>% 
  ggplot() +
  geom_hline(yintercept=1,linetype=2,alpha=.5) +
  geom_ribbon(aes(x=day,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=ara_name),alpha=.5) +
  geom_line(aes(x=day,y=exp(mean),colour=ara_name)) +
  facet_wrap(~ara_name) +
  scale_colour_discrete(guide="none") +
  scale_fill_discrete(guide="none") +
  scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
  coord_cartesian(ylim=c(.05,20)) +
  labs(title="Deviations from average time trend by ARA",x="Day",y="Relative viral load by ARA") 
#'  
#' Results are very similar to model A5.3. We find that the spatial structure explains between 0 and 70% of the spatial heterogeneity.
#' 
