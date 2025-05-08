#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: correlation regional trends hospit
# creation: jriou 
# init date: 2024-06-20
#:::::::::::::::::::::::::::::

mw_142_regional_correlation_lag = function(dat,mod,corr,lags=seq(-5,5,by=1)*7,type="hospitalizations",pprint=FALSE) {
  
  
  # Extract surveillance data -----------------------------------------------
  
  oblig = suppressMessages(readr::read_csv("../data/foph_oblig/data.csv")) 
  
  report = oblig %>% 
    dplyr::filter(!is.na(value)) %>% 
    dplyr::filter(valueCategory==type,georegion_type=="canton",agegroup=="all",sex=="all",testResult=="positive") %>% 
    dplyr::mutate(date=ISOweek::ISOweek2date(paste0(temporal,"-4"))) %>% 
    dplyr::left_join(dplyr::select(corr,georegion=kt,NUTS2_name) %>%  distinct(),by = join_by(georegion)) %>% 
    dplyr::group_by(NUTS2_name,date) %>% 
    dplyr::summarise(count=sum(value),pop=sum(pop),.groups="drop") %>%
    dplyr::mutate(inc=count/pop*1e5) %>% 
    dplyr::filter(!is.na(inc))
  
  if(FALSE) {
    ggplot(report) +
      geom_line(aes(x=date,y=inc,colour=NUTS2_name))
  }
  
  
  # Extract model time trend ------------------------------------------------
  
  alldays = unique(dat$day)
  alldays = alldays[order(alldays)]
  ndays = length(alldays)
  corr_periods = dat %>% 
    dplyr::select(date,period) %>% 
    dplyr::distinct()
  corr_days = tibble(date=seq.Date(from=min(dat$date),to=max(dat$date)-1,by=1)) %>% 
    mutate(day=row_number()-1) %>% 
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
    filter(!is.na(mean),!is.na(date)) %>% 
    left_join(report,by = join_by(NUTS2_name, date))
  
  if(FALSE) {
    ggplot(tt,aes(x=date)) +
      geom_ribbon(aes(ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
      geom_line(aes(y=exp(mean),colour=NUTS2_name),linewidth=.9) +
      geom_point(aes(y=inc)) +
      facet_wrap(~NUTS2_name) +
      coord_cartesian(ylim=c(0,20))
  }
  
  allcor_lag = NULL
  for(ll in lags) {
    if(ll<0) {
      tmp = tt %>% 
        group_by(NUTS2_name,period) %>% 
        mutate(lag_inc=lag(inc,abs(ll))) %>% 
        filter(!is.na(lag_inc))
    } else if (ll>=0) {
      tmp = tt %>% 
        group_by(NUTS2_name,period) %>% 
        mutate(lag_inc=lead(inc,abs(ll))) %>% 
        filter(!is.na(lag_inc))
    }
    tmp = tmp %>% 
      summarise(cor=cor(lag_inc,mean,method="pearson"),
                .groups="drop") %>% 
      mutate(lag=ll,
             colcor=ifelse(cor>.7,"white","black"))
    allcor_lag = bind_rows(allcor_lag,tmp)
  }
  
  g= allcor_lag %>%
    ggplot() +
    geom_line(aes(x=lag,y=cor,colour=NUTS2_name)) +
    facet_wrap(~period) +
    geom_vline(xintercept=0,linetype=2) +
    labs(x="Delay for hospitalizations",y="Correlation",title=NULL, colour="Region") +
    theme(legend.position="right")
  
  sum_allcor = allcor_lag %>% 
    group_by(period,NUTS2_name) %>% 
    filter(cor==max(cor),period %in% c(1,2,3)) %>% print(n=50)
  sum_allcor
  return(g)
  
}