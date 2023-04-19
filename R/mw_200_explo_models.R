#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: exploratory models
# creation: jriou
# init date: 2023-03-08
#:::::::::::::::::::::::::::::

mw_200_explo_models = function(ww) {
  
  # A - simple regression on one ARA ----------------------------------------
  
  # select one ARA
  chosen_ara = "Aarau"
  ww_one = ww %>% 
    filter(ara_name==chosen_ara) %>% 
    # remove missing
    filter(!is.na(vl)) %>%
    # log
    mutate(logvl=log(vl),
           logrep=log(reported_cases+1)) %>% 
    # manage dates
    mutate(day=as.numeric(date)-min(as.numeric(date))) %>% 
    rownames_to_column()
  
  # plots
  ppp_day = function(dat,mod) {
    mod$summary.fitted.values %>% 
      dplyr::bind_cols(dat) %>% 
      ggplot(aes(x=date)) +
      geom_point(aes(y=logvl)) +
      geom_ribbon(aes(ymin=`0.025quant`,ymax=`0.975quant`),alpha=.5,fill=cust_cols[1]) +
      geom_line(aes(y=mean),colour=cust_cols[1]) +
      annotate(x=min(dat$date),y=min(dat$logvl),geom="text",label=paste0("WAIC=",round(mod$waic$waic)),hjust=0)
  }
  
  # model A1: gamma regression
  ma1 = INLA::inla(logvl ~ 1 + day, 
                   data = ww_one, 
                   family = "gamma",
                   control.compute = list(waic=TRUE))
  summary(ma1)
  ppp_day(ww_one,ma1)
  ma1$waic$waic
  
  # model A2: splines
  dknots = seq(from=0,to=max(ww_one$day),by=30)
  ma2 = INLA::inla(logvl ~ 1 + 
                     bs(day,knots=dknots), 
                   data = ww_one, 
                   family = "gamma",
                   control.compute = list(waic=TRUE))
  summary(ma2)
  ppp_day(ww_one,ma2)
  ma2$waic$waic
  
  # model A3: rw1
  ma3 = INLA::inla(logvl ~ 1 + 
                     f(day,
                       model="rw1",
                       hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))), 
                   data = ww_one, 
                   family = "gamma",
                   control.compute = list(waic=TRUE))
  summary(ma3)
  ppp_day(ww_one,ma3)
  
  # model A4: rw2
  ma4 = INLA::inla(logvl ~ 1 + 
                     f(day,
                       model="rw2",
                       hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))), 
                   data = ww_one, 
                   family = "gamma",
                   control.compute = list(waic=TRUE))
  summary(ma4)
  ppp_day(ww_one,ma4)
  
  # model A5: rw2 + method change
  ma5 = INLA::inla(logvl ~ 1 + 
                     f(day,
                       model="rw2",
                       hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                     method, 
                   data = ww_one, 
                   family = "gamma",
                   control.compute = list(waic=TRUE))
  summary(ma5)
  ppp_day(ww_one,ma5)
  
  # model A6: reported cases
  cowplot::plot_grid(ggplot(ww_one) + geom_point(aes(x=date,y=logrep)),
                     ggplot(ww_one) + geom_point(aes(x=date,y=logvl),col="red"),
                     ncol=1)
  ma6 = INLA::inla(logvl ~ 1 + logrep, 
                   data = ww_one, 
                   family = "gamma",
                   control.compute = list(waic=TRUE))
  summary(ma6)
  ppp_day(ww_one,ma6)
  
  # model A7: reported cases by period (defined using low points in reported cases)
  cut_date = ymd("2022-05-15","2022-09-15","2023-01-01")
  cut_day = as.numeric(cut_date) - min(as.numeric(ww_one$date))
  ww_one = ww_one %>% 
    dplyr::mutate(period=cut(day,breaks=c(-1,cut_day,max(ww_one$day)+1)))
  cowplot::plot_grid(ggplot(ww_one) + geom_point(aes(x=date,y=logrep),col="red") + geom_vline(xintercept=cut_date,linetype=2),
                     ggplot(ww_one) + geom_point(aes(x=date,y=logvl)) + geom_vline(xintercept=cut_date,linetype=2),
                     ncol=1)
  ggplot(ww_one) + geom_point(aes(x=logrep,y=logvl,col=period)) + facet_wrap(~period)
  ma7 = INLA::inla(logvl ~ 1 + 
                     logrep:period, 
                   data = ww_one, 
                   family = "gamma",
                   control.compute = list(waic=TRUE))
  summary(ma7)
  ppp_day(ww_one,ma7) + geom_vline(xintercept=cut_date,linetype=2)
  
  ww_one %>% 
    bind_cols(ma7$summary.fitted.values) %>% 
    ggplot(aes(x=logrep)) + 
    geom_point(aes(y=logvl,col=period)) + 
    geom_ribbon(aes(ymin=`0.025quant`,ymax=`0.975quant`),alpha=.5,fill=cust_cols[1]) +
    geom_line(aes(y=mean),colour=cust_cols[1]) +
    facet_wrap(~period)
  
  
  
  ww_one %>% 
    mutate(mean=predict(ma7b)) %>% 
    ggplot(aes(x=logrep)) + 
    geom_point(aes(y=logvl,col=period)) + 
    # geom_ribbon(aes(ymin=`0.025quant`,ymax=`0.975quant`),alpha=.5,fill=cust_cols[1]) +
    geom_line(aes(y=mean),colour=cust_cols[1]) +
    facet_wrap(~period)
  
  
  # model A8: reported cases by period + remaining time trend
  ma8 = INLA::inla(logvl ~ 1 + 
                     logrep:period +
                     f(day,period,
                       model="rw2",
                       hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.001)))), 
                   data = ww_one, 
                   family = "gaussian",
                   control.compute = list(waic=TRUE))
  summary(ma8)
  ppp_day(ww_one,ma8)
  
  
}
