#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: correlation regional trends hospit
# creation: jriou 
# init date: 2024-06-20
#:::::::::::::::::::::::::::::

mw_140_regional_correlation = function(dat,mod,corr,type="hospitalizations",pprint=FALSE) {
  

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
      facet_wrap(~NUTS2_name)
    tt %>% 
      # filter(NUTS2_name=="Zurich") %>% 
      ggplot(aes(x=inc)) +
      geom_point(aes(y=exp(mean),colour=period)) +
      # geom_smooth(method='lm') +
      facet_grid(NUTS2_name~period,scales="free") +
      scale_y_log10()
  }
  allcor = tt %>% 
    filter(!is.na(inc)) %>% 
    group_by(NUTS2_name,period) %>% 
    summarise(cor=cor(inc,mean,method="pearson"),.groups="drop") %>% 
    mutate(colcor=ifelse(cor>.7,"white","black"))
  
  g = allcor %>% 
    ggplot() +
    geom_tile(aes(x=period,y=NUTS2_name,fill=cor)) +
    geom_text(aes(x=period,y=NUTS2_name,label=sprintf("%.2f",cor),colour=colcor),size=3.4) +
    scale_fill_distiller(palette = "Greys",direction = 1,limits=c(0,1)) +
    scale_y_discrete(expand=c(0,0),limits=rev) +
    scale_x_discrete(expand=c(0,0)) +
    scale_colour_manual(values=c("black","white"),guide="none") +
    labs(x="Phase",y="Region",fill="Correlation    ") +
    theme(legend.position="bottom")
  
  if(pprint) {
    g=allcor
  }
  if(FALSE) {
    tt %>% 
      filter(!is.na(inc)) %>% 
      ggplot(aes(x=inc,y=exp(mean))) +
      geom_point() +
      geom_smooth(method='lm',formula = 'y ~ x',aes(colour=NUTS2_name)) +
      facet_grid(NUTS2_name~period) +
      scale_colour_discrete(guide="none") +
      scale_y_log10() +
      geom_text(data=allcor,aes(label=paste(sprintf("%.2f",cor)," "),x=Inf, y=0.25,
                                vjust=1,hjust=1)) +
      labs()
  }
  return(g)
  
}