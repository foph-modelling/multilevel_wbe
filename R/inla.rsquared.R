#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: helper functions to manage inla output in model_dev_A
# creation: jriou 
# init date: 2024-06-19
#:::::::::::::::::::::::::::::


inla.rsquared = function(dat,mod) {
  lims = dat %>% 
    dplyr::filter(below_lod==0,below_loq==0) %>% 
    dplyr::summarise(minvl=min(vl),maxvl=max(vl))
  tt = mod$summary.fitted.values %>% 
    dplyr::bind_cols(dat)
  r2 = 1 - sum((tt$vl - tt$mean)^2) / sum((tt$vl - mean(tt$vl))^2)
  return(r2)
} 