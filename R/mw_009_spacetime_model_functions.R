library(sf)
library(INLA)
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
library(utils)

fit_inla_model = function(wwdata, 
                          catchment_centroids, 
                          plz_pos,
                          global_temp_model='rw1', 
                          sptm_temp_model='ar1', 
                          fam='gamma',
                          logvl=FALSE, 
                          save.point = NULL,
                          covariates = c('u20', 'o65', 'nec', 'pop_dens', 'lab_method')
                          )
{
  if(is.null(save.point)){
    save.point = paste0('outputs/last_run_', Sys.time())
    dir.create(save.point)
  }
  small = 1e-23 # small value to offset 0s
  
  wwdata = wwdata[order(day, ara_id),]
  # Merge catchment centroids to WW measurements
  catchment_centroids_ww = catchment_centroids %>% arrange(day, ara_id)

  #generate data input for INLA model
  cnetroid_coords = st_coordinates(catchment_centroids_ww)
  d = cbind(data.table(cnetroid_coords), data.table(st_drop_geometry(catchment_centroids_ww))[, c('ara_id', 'day', 'pop_total', 'vl_stand', ..covariates)])
  
  #generate container for predicted values of INLA model - used after run (to be moved to own function )
  #plzcoords = st_coordinates(plz_pos)
  #pcoords = data.table()
  #for(day in unique(d$day)){
  #  pcoords_date = cbind(plzcoords, day, plz_pos$PLZ)
  #  pcoords = rbind(pcoords, pcoords_date)
  #}
  #names(pcoords) <- c("x", "y", "time", "PLZ")
  #pcoords = merge(unique(pcoords), unique(covariate_pops_plz_norm), by=c('PLZ'), how='left')
  
  
  # generate mesh for INLA model 
  max.edge = diff(range(cnetroid_coords[,1]))/(3*5)
  mesh <- inla.mesh.2d(
    loc = d[,c('X', 'Y')], max.edge = max.edge * c(1,2), cutoff=1000
  )
  plot(mesh)
  points(cnetroid_coords, col = "red")
  
  # generate stochastic pdes for inla model
  spde <- inla.spde2.pcmatern(
    mesh = mesh, alpha = 2,
    prior.range = c(100000, 0.01), # P(range < 10000) = 0.01
    prior.sigma = c(0.5, 0.01) # P(sigma > 0.5) = 0.01
  )
  
  
  # set indices for INLA model 
  group <- d$day - min(d$day, na.rm = T) + 1
  timesn <- length(unique(group))
  indexs <- inla.spde.make.index("s",
                                 n.spde = spde$n.spde,
                                 n.group = timesn)
  
  # generate A and Ap for INLA model 
  A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(d[,c('X', 'Y')]), group = group)
  #Ap <- inla.spde.make.A(mesh = mesh, loc = as.matrix(pcoords[,c(2, 3)], ncol=2), group=as.matrix(pcoords[,4] ))
  
  
  
  if(logvl==TRUE){
    d$vl_stand = log(d$vl_stand + small)
  }
  else{
    d$vl_stand = d$vl_stand + small
  }
  
  # construct INLA stacks
  stk.e <- inla.stack(
    tag = "est",
    data = list(y = d$vl_stand),
    A = list(1, A),
    effects = list(data.frame(b0 = rep(1, nrow(d)), u20=d$u20, o65=d$o65, nec=d$nec, pop_dens=d$pop_dens, time.index=group, lab_method = d$lab_method), s = indexs)
  )
  
  
  # set priors and link covariate matrix (set to 1)
  
  #rprior3 <- list(theta = list(prior = "pc.cor1", param = c(0.9, 0.9)))
  tprior_global <- 'list(prec = list(prior = "pc.prec", param = c(1, 0.01)))'
  tprior_local <- 'list(theta1 = list(prior="pc.prec", param=c(1, 0.01)), rho = list(prior="pc.cor1", param=c(0.1, 0.9)))' 
  
  
  # model formula 
  
  
  #formula <- make_sptm_formula(sptm_temp_model = 'ar1', 
  #                             global_temp_model = 'rw1', 
  #                             covariates = covariates, 
  #                             tprior_global = tprior_global, 
  #                             tprior_local = tprior_local)
    #rprior3 <- list(theta = list(prior = "pc.cor1", param = c(0.9, 0.9)))
  rprior1 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  rprior3 <- list(theta1 = list(prior="pc.prec", param=c(1, 0.01)), theta = list(prior="pc.cor1", param=c(0.1, 0.9))) 
  
  
  # model formula 
  
  
  formula <- y ~ 0 + 
                 b0 + u20 + o65 + nec  + pop_dens + lab_method + 
                 f(time.index, model='rw1', hyper=rprior1) +
                 f(s, model = spde, group = s.group, control.group = list(model = "ar1", hyper = rprior3))
  # run the model 
  
  res <- inla(formula,
              family=fam,
              data = inla.stack.data(stk.e),
              control.compute = list(config=TRUE),
              control.predictor = list(
                compute = TRUE,
                link = 1,
                A = inla.stack.A(stk.e)
              ),
              safe=TRUE,
              num.threads = 8
  )
  
  
  
  marginals =data.table() 
  for(var in res$names.fixed){
    var_frame = data.table(res$marginals.fixed[[var]])
    var_frame[, var:=var]
    marginals = rbind(marginals, var_frame)
  }
  for(var in names(res$marginals.hyperpar)){
    var_frame = data.table(res$marginals.hyperpar[[var]])
    var_frame[, var:=var]
    marginals = rbind(marginals, var_frame)
  }
  
  saveRDS(list(marginals=marginals, A.indexs.spde=list(A=A,indices=indexs,spde=spde), stk=stk.e, pcoords=NULL), paste0(save.point, '/model_inputs.rds'))
  
  
  return(list(marginals=marginals, res=res, A.indexs.spde=list(A=A,indices=indexs,spde=spde), stk=stk.e, pcoords=NULL))
}



make_sptm_formula = function(global_temp_model='rw2', sptm_temp_model='ar1', covariates=c(), tprior_global, tprior_local){
  
  
  
  fixed_effects <- ''
  if (!is.null(covariates)) {
    fixed_effects <- paste0(paste(covariates,collapse=' + '),' + ')
  }
  
  f <- paste('y ~ 0 + b0 + ',                                           
             fixed_effects,                                                        
             paste0('   f(time.index,model=', global_temp_model, ', hyper =', tprior_global),  
             paste0(' + f(s, model = A.indexs.spde$spde, group = s.group,
                   control.group = list(model =', sptm_temp_model, ', hyper =',tprior_local,'))')
  )
  return(f)
}
  



specify_space_time_model <- function(fixed_effects = NULL, region_effect){
  #
  rgn_effect <- ''
  if (region_effect=='re') {   # region random effect
    rgn_effect <- ' + f(region.index,model=\'iid\',hyper = hyper.prec.region)'
  }
  fx <- ''
  if (!is.null(fixed_effects)) {
    fx <- paste(fixed_effects,collapse=' + ')
  }
  # fixed effects are a list
  f <- paste('y ~ 0 + b0 + ',                                           # intercept
             fx,                                                        # covariates/fixed effects
             ' + f(time.index,model=\'rw1\',hyper = hyper.prec.time)',  # overall time profile
             ' + f(site.index,model=\'iid\',hyper = hyper.prec.space)', # overall site-effects
             rgn_effect,                                                # region random effect
             # spatially-temporally correlated space-time interactions
             ' + f(s, model = A.indexs.spde$spde,
                   group = s.group,
                   control.group = list(model = "ar1", hyper = pcprior.rho()))'
  )
  return(f)
}


get_samples_from_inla_model = function(inla_results, covariates, pred_coords_covars, covariate_matrix, pred_times=NULL, nsims=500, model=NULL, model_dir='', suffix='', log_res=FALSE){
  
  if(suffix != ''){
    suffix = paste0('_', suffix)
  }
  
  message('constructing prediction containers')
  
  fixed_pars_cols = covariates
  fixed_pars = covariates
  
  res = inla_results[['res']]
  A.indexs.spde = inla_results[['A.indexs.spde']]
  stk = inla_results[['stk']]
  
  pcoords = pred_coords_covars

  
  locations = unique(pcoords[order(time, PLZ),]$PLZ)
  
  message(paste0('Drawing ', nsims, ' samples for predictions at ', length(time_steps), ' time points and ', length(locations), ' locations...'))
  

  
  fixed.effects.in.order = pcoords[order(time, PLZ),..fixed_pars_cols]
  
  ntimes = length(time_steps)
  
  coordinates = unique(pcoords[,c('x','y')])

  message('Sampling posterior...')

  samples = inla.posterior.sample(n=nsims, result = res)

  message('Posterior sampling complete!! \n\n')
  
  
  message('extracting intercept ...')
  intercept <- INLA::inla.posterior.sample.eval('b0',samples)[1,]
  
  message('extracting fixed effect covefficients ...')
  fixed <- array(0,c(nsims,length(fixed_pars)))
  for (i in 1:length(fixed_pars)) {
    fixed[,i] <- INLA::inla.posterior.sample.eval(fixed_pars[i],samples)[1,]
  }
  
  message('extracting hyper parameters ...')
  extract.joint.hyperpar <- function(samples,para) {
    out <- unlist(lapply(samples,function(x){x$hyperpar[para]}))
    return(out)
  }
  
  
  sigma.residual <- sqrt(1/extract.joint.hyperpar(samples,'Precision parameter for the Gamma observations'))
  rho <- extract.joint.hyperpar(samples,'GroupRho for s')
  

  
  
  sigma.matern <- extract.joint.hyperpar(samples,'Stdev for s')
  
  range <- extract.joint.hyperpar(samples,'Range for s')
  
  message('extracting global temporal parameters ... ')
  time <- t(INLA::inla.posterior.sample.eval('time.index',samples))
  pred_time_index <- rep(1:ntimes,each=nrow(pcoords[time==1,]))
  
  message('extracting spatio-temporal parameters ...')
  cont <- attr(samples, which = ".contents", exact = TRUE)
  id <- which(cont$tag=='s')  #  's' corresponds to the name associated with the space-time interaction in the fitted model: f(s, model = A.indexs.spde$spde ...
  start <- cont$start[id]
  end <- cont$length[id] + start - 1
  ss <- lapply(samples,function(x){x$latent[start:end,1]})
  
  ss.array <- t(matrix(unlist(ss), ncol=nsims))
  
  #extract.1d.sptm.interactions <- function(interactions,itm=NULL,isp=NULL,site.index,time.index) {
  #  if (!is.null(itm) & !is.null(isp)) {
  #    stop('ERROR in extract.1d.sptm.interactions: This function can only extract one dimension, space or time but not both')
  #  }
  #  if (!is.null(itm)) {
  #    # extracting all terms at time itm
  #    ids <- which(time.index==itm)
  #    if (is_consecutive(site.index[ids],unique=FALSE)) {
  #      out <- interactions[,ids]
  #    } else {
  #      stop('ERROR in extract.1d.sptm.interactions: order of site indices is wrong')
  #    }
  #  }
  #  if (!is.null(isp)) {
  #    # extracting all terms associated with site isp
  #    ids <- which(site.index==isp)
  #    if (is_consecutive(time.index[ids],unique=FALSE)) {
  #      out <- interactions[,ids]
  #    } else {
  #      stop('ERROR in extract.1d.sptm.interactions: order of time indices is wrong')
  #    }
  #  }
  #  return(out)
  #}
  #
  #ntimes <- 10 #ncol(sims.list.joint$time)
  #s.by.time <- lapply(1:ntimes,function(it) {
  #  extract.1d.sptm.interactions(ss.array,itm=it,isp=NULL,
  #                               site.index=A.indexs.spde$indices$s,
  #                               time.index=A.indexs.spde$indices$s.group)
  #})
  
  s.by.time <- lapply(1:ntimes,function(it) {
    ids = which(A.indexs.spde$indices$s.group==it)
    ss.array[,ids]}
  )
  
  
  l = length(s.by.time)
  d = dim(s.by.time[[1]])
  s.by.time.a <- array(0,c(l,d[1],d[2]))
  for (i in 1:l) {
    s.by.time.a[i,,] <- s.by.time[[i]]
  }
  
  s.by.time = s.by.time.a
  
  pred.mesh <- INLA::inla.mesh.projector(A.indexs.spde$spde$mesh,loc=as.matrix(coordinates))
  
  eta.plz.sims <- lapply(1:nsims,function(isim) {
    x <- s.by.time[,isim,]
    #  this projects onto a (nlsoas x ntimes) matrix of eta values over the LSOA centroids
    m <- apply(x,1,function(vs){INLA::inla.mesh.project(pred.mesh,vs)})
    #  format the matrix to long in the form of eta_11,..., eta_N1,eta_12,....,eta_N2,...,eta_1T,..., eta_NT
    #  eg c(matrix(1:12,nrow=3,byrow=FALSE))
    c(m)
  })
  
  message('extracting error ...')
  
  eta <- t(matrix(unlist(eta.plz.sims), ncol=nsims))
  
  eps <- sapply(1:nsims,function(x){rnorm(nrow(pcoords),0,sigma.residual[x])})  # per space-time cell
  epsilon  = t(eps)
  
  
  message('calculating predicted values ... ')
  
  sims.pred <- intercept +
    fixed %*% t(fixed.effects.in.order)+ 
    time[,pred_time_index] +
    #vsp +
    eta +
    epsilon
  
  
  if(log_res==TRUE){
    sims.pred.t = exp(t(sims.pred))
  }else{
    sims.pred.t = t(sims.pred)
  }
  
  message(paste0('\n \nsaving predictions at: ', model_dir,'/posterior_predictions_', as.character(nsims), '_', model, suffix, '.rds'))
  
  saveRDS(object = sims.pred.t, file=paste0(model_dir,'/posterior_predictions_', as.character(nsims), '_', model, suffix, '.rds' ))
  saveRDS(object = list(intercept, fixed.effects.in.order, fixed, eta, epsilon),file=paste0(model_dir,'/posterior_params_', as.character(nsims), '_', model, suffix, '.rds' ))
  
  
  
}


library(sf)


get_plz_areas = function(file_path = 'data/spatial/plz/PLZO_SHP_LV95/PLZO_PLZ.shp', crs_required=NULL) {
  # read file
  plz_area = st_read(file_path)
  
  # set to required crs
  if(!is.null(crs_required)){
    plz_area = st_transform(plz_area, crs = crs_required)
  }
  
  #merge split polygons
  plz_area = rmapshaper::ms_dissolve(plz_area,'PLZ')
  
  
  #return shapes in sf object
  plz_area
  
}

get_plz_centroids = function(file_path = 'data/spatial/plz/PLZO_SHP_LV95/PLZO_PLZ.shp', crs_required=NULL){
  
  #load plz areas
  plz_area = get_plz_areas(file_path, crs_required)
  
  #calculate centroids
  plz_pos = st_centroid(plz_area)
  
  #return centroids in sf object
  plz_pos
}





