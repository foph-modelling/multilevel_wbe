#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: helper functions to manage inla output in model_dev_A
# creation: jriou 
# init date: 2023-05-09
#:::::::::::::::::::::::::::::


# k-fold cross validation ----
kfoldcv_vl = function(dat,mod,k=10,showplot=FALSE) {
  # extract formula
  ff = as.character(mod$.args$formula)
  ff = as.formula(paste(ff[2],ff[1],ff[3]))
  # adapt data
  ser = rep(1:k,length.out=nrow(dat))
  ser = ser[order(ser)]
  dd = dat %>% 
    # dplyr::mutate(k=sample(1:k,nrow(dat),replace=TRUE),
    dplyr::mutate(k=ser,
                  mean=NA,lwb=NA,upb=NA)
  # run k-fold cv
  for(i in 1:k) {
    rr = which(dd$k==i)
    ddx = dd %>% 
      dplyr::mutate(vl = ifelse(k==i,NA,vl)) 
    mm = INLA::inla(ff,
                    data = ddx,
                    family = "gamma",
                    control.mode = list(result = mod, restart=TRUE),
                    control.predictor = list(compute=TRUE,link=1))
    dd[rr,"mean"] = mm$summary.fitted.values[rr,"mean"]
    dd[rr,"lwb"] = mm$summary.fitted.values[rr,"0.025quant"]
    dd[rr,"upb"] = mm$summary.fitted.values[rr,"0.975quant"]
  }
  if(showplot==TRUE) {
    g = dd %>% 
      filter(!is.na(vl)) %>% 
      ggplot(aes(x=date)) +
      geom_point(aes(y=vl)) +
      geom_ribbon(aes(ymin=lwb,ymax=upb),alpha=.5,fill=cust_cols[2]) +
      geom_line(aes(y=mean),colour=cust_cols[2]) +
      scale_y_continuous(trans="log",
                         breaks = trans_breaks("log", function(x) exp(x)),
                         labels = trans_format("log", math_format(e^.x))) 
    print(g)
  }
  oo = dd %>% 
    dplyr::filter(!is.na(vl)) %>% 
    dplyr::transmute(se=(vl-mean)^2,
                     coverage=ifelse(vl>lwb & vl<upb,1,0),
                     logse=(log(vl)-log(mean))^2,
                     logsharpness=(log(upb)-log(lwb))/2) %>% 
    summarise(rmse=sqrt(mean(se)),
              coverage=mean(coverage),
              logrmse=sqrt(mean(logse)),
              logsharpness=mean(logsharpness)) %>% 
    unlist() %>% unname()
  output = list(RMSE=round(oo[1],2),
                coverage=round(oo[2],2),
                logRMSE=round(oo[3],2),
                sharpness=round(oo[4],2))
  return(output)
}

avg_time_trend = function(dat,mod) {
  lims = dat %>% 
    dplyr::filter(below_lod==0,below_loq==0) %>% 
    dplyr::summarise(minvl=min(vl),maxvl=max(vl))
  alldays = 0:(max(dat$day)-1)
  ndays = length(alldays)
  corr_periods = dat %>% 
    dplyr::select(date,period) %>% 
    dplyr::distinct()
  corr_days = tibble(day=alldays,
                     date=seq.Date(from=min(dat$date),to=max(dat$date)-1,by=1)) %>% 
    left_join(corr_periods,by = join_by(date))
  tt = mod$summary.random$day %>% 
    dplyr::mutate(day=ID) %>% 
    dplyr::left_join(corr_days,by = join_by(day)) %>% 
    as_tibble() %>% 
    filter(!is.na(mean),!is.na(date))
  g = ggplot(tt,aes(x=date)) +
    geom_hline(yintercept=1,linetype=2,alpha=.5) +
    geom_ribbon(aes(ymin=exp(`0.025quant`),ymax=exp(`0.975quant`)),alpha=.5) +
    geom_line(aes(y=exp(mean))) +
    scale_y_continuous(trans="log",breaks=c(0,.2,.5,1,2,5,10)) +
    labs(x="Date",y="Relative VL",title="Average time trend")
  return(g)
}


avg_time_trend_reg = function(dat,mod) {
  lims = dat %>% 
    dplyr::filter(below_lod==0,below_loq==0) %>% 
    dplyr::summarise(minvl=min(vl),maxvl=max(vl))
  alldays = 0:(max(dat$day)-1)
  ndays = length(alldays)
  corr_periods = dat %>% 
    dplyr::select(date,period) %>% 
    dplyr::distinct()
  corr_days = tibble(day=alldays,
                     date=seq.Date(from=min(dat$date),to=max(dat$date)-1,by=1)) %>% 
    left_join(corr_periods,by = join_by(date))
  allnuts2 = dat %>% 
    dplyr::select(NUTS2,NUTS2_name) %>% 
    dplyr::distinct() %>% 
    dplyr::arrange(NUTS2)
  tt = mod$summary.random$day %>% 
    dplyr::mutate(day=rep(alldays,nrow(allnuts2)),
                  NUTS2=rep(allnuts2$NUTS2,each=ndays)) %>% 
    dplyr::left_join(allnuts2,by = join_by(NUTS2)) %>% 
    dplyr::left_join(corr_days,by = join_by(day)) %>% 
    as_tibble() %>% 
    filter(!is.na(mean),!is.na(date))
  g = ggplot(tt,aes(x=date)) +
    geom_hline(yintercept=1,linetype=2,alpha=.5) +
    # geom_ribbon(aes(ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
    geom_line(aes(y=exp(mean),colour=NUTS2_name)) +
    scale_y_continuous(trans="log",breaks=c(0,.2,.5,1,2,5,10)) +
    labs(x="Date",y="Relative VL",title="Average time trend by region")
  return(g)
}

# plot posterior predictive plot to viral load data ----
ppp_vl = function(dat,mod) {
  lims = dat %>% 
    dplyr::filter(below_lod==0,below_loq==0) %>% 
    dplyr::summarise(minvl=min(vl),maxvl=max(vl))
  mod$summary.fitted.values %>% 
    dplyr::bind_cols(dat) %>% 
    ggplot(aes(x=date)) +
    geom_point(aes(y=vl)) +
    geom_ribbon(aes(ymin=`0.025quant`,ymax=`0.975quant`),alpha=.5,fill=cust_cols[1]) +
    geom_line(aes(y=mean),colour=cust_cols[1]) +
    scale_y_continuous(trans="log",
                       breaks = trans_breaks("log", function(x) exp(x)),
                       labels = trans_format("log", math_format(e^.x))) +
    annotate(x=min(dat$date),y=lims$minvl,geom="text",label=paste0("WAIC=",round(mod$waic$waic)),hjust=0) +
    coord_cartesian(ylim=c(lims$minvl,lims$maxvl)) +
    labs(x="Date",y="Viral load",title="Posterior predictive plot")
}
ppp_vl_ara = function(dat,mod) {
  lims = dat %>% 
    dplyr::filter(below_lod==0,below_loq==0) %>% 
    dplyr::summarise(minvl=min(vl),maxvl=max(vl))
  mod$summary.fitted.values %>% 
    dplyr::bind_cols(dat) %>% 
    ggplot(aes(x=date)) +
    geom_point(aes(y=vl),alpha=.3) +
    geom_ribbon(aes(ymin=`0.025quant`,ymax=`0.975quant`),alpha=.5,fill=cust_cols[1]) +
    geom_line(aes(y=mean),colour=cust_cols[1]) +
    scale_y_continuous(trans="log",
                       breaks = trans_breaks("log", function(x) exp(x)),
                       labels = trans_format("log", math_format(e^.x))) +
    coord_cartesian(ylim=c(lims$minvl,lims$maxvl)) +
    labs(x="Date",y="Viral load",title="Posterior predictive plot") + 
    facet_wrap(~ara_name) + 
    theme(axis.text.x = element_text(angle=45,hjust=1))
}

# extract and exponent model parameters ----
summary_exp_vl = function(mod, pars, order=FALSE) {
  o = mod$summary.fixed %>% 
    rownames_to_column() %>% 
    as_tibble() %>% 
    dplyr::filter(grepl(pars,rowname)) %>% 
    dplyr::transmute(rowname=rowname,
                     `exp(beta)`=round(exp(mean),2),
                     `0.025quant`=round(exp(`0.025quant`),2),
                     `0.975quant`=format(round(exp(`0.975quant`),2),scientific=FALSE)) %>% 
    column_to_rownames()
  if(order) {
    o = o %>% 
      dplyr::arrange(-`exp(beta)`)
  }
  return(o)
}