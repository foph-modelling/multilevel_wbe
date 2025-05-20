#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: hospits by period by region
# creation: jriou 
# init date: 2025-05-08
#:::::::::::::::::::::::::::::

mw_143_regional_hospits = function(dat,type="hospitalizations") {
  

# Extract surveillance data -----------------------------------------------

  oblig = suppressMessages(readr::read_csv("../data/foph_oblig/data.csv")) 
  cut_dates = lubridate::ymd(c(min(dat$date)-1,controls$period_dates),max(dat$date)+1)
  
  report = oblig %>% 
    dplyr::filter(!is.na(value)) %>% 
    dplyr::filter(valueCategory==type,georegion_type=="canton",agegroup=="all",sex=="all",testResult=="positive") %>% 
    dplyr::mutate(date=ISOweek::ISOweek2date(paste0(temporal,"-4"))) %>% 
    dplyr::left_join(dplyr::select(corr,georegion=kt,NUTS2_name) %>%  distinct(),by = join_by(georegion)) %>% 
    dplyr::group_by(NUTS2_name,date) %>% 
    dplyr::summarise(count=sum(value),pop=sum(pop),.groups="drop") %>%
    dplyr::mutate(inc=count/pop*1e5) %>% 
    dplyr::filter(!is.na(inc)) %>% 
    dplyr::mutate(period=cut(date,breaks=cut_dates),
                  period=as.factor(as.numeric(as.factor(period)))) %>% 
    dplyr::filter(!is.na(period))
    
  g = report %>% 
    dplyr::group_by(period,NUTS2_name) %>% 
    dplyr::summarise(cuminc=sum(inc)) %>% 
    ggplot() +
    geom_col(aes(x=period,y=cuminc,fill=NUTS2_name),colour="black",position = position_dodge()) +
    labs(x="Period",y="Hospitalisation rate (per 100,000)",fill="Region")
  
  
  return(g)
  
}