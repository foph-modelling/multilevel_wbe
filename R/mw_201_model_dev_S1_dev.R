controls = readRDS(fs::path("../savepoints/savepoint_2023-05-15/controls.rds"))

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

ww_all = ww_all %>% filter(day1<60 & NUTS2_name=='Mittelland')
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


mw_100_desc_table(ww_all) %>%
  dplyr::mutate(across(everything(),as.character)) %>%
  tidyr::gather() %>%
  flextable::flextable(cwidth=c(4,4))

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

#ww_all = ww_all[order(day, ara_id), ]

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
                              hol=ww_all$hol, 
                              ara1=ww_all$ara1, 
                              day1=ww_all$day1, 
                              ara2=ww_all$ara2), 
                   s = indexs)
  )




#rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.9)))


# model formula with auto-regressive component - same as previous spatial model 

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



ms5.3.1 = inla(formula,
               family='gaussian',
               data = inla.stack.data(stk.e),
               control.compute = list(config=TRUE),
               control.predictor = list(
                 compute = TRUE,
                 A = inla.stack.A(stk.e)
               ), 
               safe=TRUE
                 
  )


dp_s = copy(ww_all)

index_s <- inla.stack.index(stack = stk.e, tag = "est")$data
dp_s$pred_mean = ms5.3.1$summary.fitted.values[index_s, "mean"]
dp_s$pred_lower = ms5.3.1$summary.fitted.values[index_s, "0.025quant"]
dp_s$pred_upper = ms5.3.1$summary.fitted.values[index_s, "0.975quant"]

ggplot(dp_s[below_lod==0 & below_loq==0, ]) + 
  geom_point(aes(x=day1, y=vl_stand), size=0.3) + 
  geom_line(aes(x=day1, y=pred_mean), color='red', linewidth=0.2, alpha=0.8)+
  geom_ribbon(aes(x=day1, ymin=pred_lower, ymax=pred_upper), fill='red', alpha=0.2)+
  facet_wrap(~ara_id)+
  scale_y_continuous(trans='log')

# model formula with random walk - the same as previous non-spatial model 

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



# run the model with Gaussian likelihood - allows negative values (consistent with Blangiardo)

ms5.3.2 = inla(formula,
               family='gaussian',
               data = inla.stack.data(stk.e),
               control.compute = list(config=TRUE),
               control.predictor = list(
                 compute = TRUE,
                 A = inla.stack.A(stk.e)
               )
               
)


dp_s = copy(ww_all)

index_s <- inla.stack.index(stack = stk.e, tag = "est")$data
dp_s$pred_mean = ms5.3.2$summary.fitted.values[index_s, "mean"]
dp_s$pred_lower = ms5.3.2$summary.fitted.values[index_s, "0.025quant"]
dp_s$pred_upper = ms5.3.2$summary.fitted.values[index_s, "0.975quant"]


ms2_plot = ggplot(dp_s[below_lod==0 & below_loq==0, ]) + 
  geom_point(aes(x=day1, y=vl_stand), size=0.3) + 
  geom_line(aes(x=day1, y=pred_mean), color='red', linewidth=0.2, alpha=0.8)+
  geom_ribbon(aes(x=day1, ymin=pred_lower, ymax=pred_upper), fill='red', alpha=0.2)+
  facet_wrap(~ara_id)+
  scale_y_continuous(trans='log')


# Fit the model with gamma likelihood - more principled for positive real values 

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


formula = y ~ 1 + beta +  
  f(below_loq,model="iid") +
  f(below_lod,model="iid") +
  f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
  f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
  f(hol,model="linear",mean.linear=0,prec.linear=.2) +
  lab_method + 
  f(s, 
    model = spde, 
    group = s.group, 
    control.group = list(model = "ar1", hyper = list(theta = list(prior = "pccor1", param = c(0, 0.9)))))



# run the model with Gaussian likelihood - allows negative values (consistent with Blangiardo)

ms5.3.5 = inla(formula,
               family='gaussian',
               data = inla.stack.data(stk.e),
               control.compute = list(config=TRUE),
               control.predictor = list(
                 compute = TRUE,
                 A = inla.stack.A(stk.e)
               )
               
)


dp_s = copy(ww_all)

index_s <- inla.stack.index(stack = stk.e, tag = "est")$data
dp_s$pred_mean = ms5.3.5$summary.fitted.values[index_s, "mean"]
dp_s$pred_lower = ms5.3.5$summary.fitted.values[index_s, "0.025quant"]
dp_s$pred_upper = ms5.3.5$summary.fitted.values[index_s, "0.975quant"]


ms5_plot = ggplot(dp_s[below_lod==0 & below_loq==0, ]) + 
  geom_point(aes(x=day1, y=vl_stand), size=0.3) + 
  geom_line(aes(x=day1, y=pred_mean), color='red', linewidth=0.2, alpha=0.8)+
  geom_ribbon(aes(x=day1, ymin=pred_lower, ymax=pred_upper), fill='red', alpha=0.2)+
  facet_wrap(~ara_id)+
  scale_y_continuous(trans='log')

