#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: correlation regional trends hospit
# creation: jriou 
# init date: 2024-06-20
#:::::::::::::::::::::::::::::

mw_141_crude_correlation = function(dat,corr,type="hospitalizations") {
  

# Extract surveillance data -----------------------------------------------

  oblig = suppressMessages(readr::read_csv("data/foph_oblig/data.csv")) 
  
  report = oblig %>% 
    dplyr::filter(!is.na(value)) %>% 
    dplyr::filter(valueCategory==type,georegion_type=="canton",agegroup=="all",sex=="all",testResult=="positive") %>% 
    dplyr::mutate(date=ISOweek::ISOweek2date(paste0(temporal,"-1"))) %>% 
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
 
  allcor = dat %>% 
    dplyr::select(wwtp_index,NUTS2_name,date,period,vl) %>% 
    dplyr::left_join(report,by = join_by(NUTS2_name, date)) %>% 
    dplyr::filter(!is.na(inc)) %>% 
    group_by(NUTS2_name,period,wwtp_index) %>% 
    summarise(cor=cor(inc,vl,method="pearson"),.groups="drop_last") %>% 
    summarise(cor=mean(cor,na.rm=TRUE),.groups="drop_last")
  
  g = allcor %>% 
    ggplot() +
    geom_tile(aes(x=period,y=NUTS2_name,fill=cor)) +
    geom_text(aes(x=period,y=NUTS2_name,label=sprintf("%.2f",cor)),size=3.4) +
    scale_fill_distiller(palette = "Greys",direction = 1.,limits=c(0,1)) +
    scale_y_discrete(expand=c(0,0),limits=rev) +
    scale_x_discrete(expand=c(0,0)) +
    labs(x="Period",y="Region",fill="Correlation") +
    theme(legend.position="bottom")
  
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