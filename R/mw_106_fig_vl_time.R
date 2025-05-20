#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for wastewater viral load
# creation: jriou 
# init date: 2023-08-07
#:::::::::::::::::::::::::::::

mw_106_fig_vl_time = function(ww) {
  plab = tibble(date=c(min(ww$date), controls$period_dates,max(ww$date))) %>% 
    mutate(date2=lead(date)) %>% 
    filter(!is.na(date2)) %>% 
    mutate(dif=(date2-date)/2) %>% 
    mutate(pos=date+dif) %>% 
    mutate(lab=as.character(row_number()),
           y=7e14)
  xlabs = c("2022-05-16","2022-09-05","2023-01-02","2023-07-03")
  g = ww %>% 
    dplyr::filter(below_lod==0,below_loq==0) %>% 
    # dplyr::group_by(ara_kt,week) %>% 
    # dplyr::summarise(vl=mean(vl,na.rm=TRUE),.groups="drop") %>% 
    ggplot() +
    geom_line(aes(x=date,y=vl,colour=vl,group=ara_id),alpha=.5) +
    scale_colour_viridis_c(trans="log",guide = "none") +
    scale_y_continuous(trans="log",breaks=c(1e10,1e11,1e12,1e13,1e14,1e15)) +
    # breaks = trans_breaks("log", function(x) exp(x)),
    # labels = trans_format("log", math_format(e^.x))) +  
    labs(x="Time",
         y="Viral load",
         fill="Mean viral load") +
    scale_x_date(expand=c(0,0),
                 breaks="4 months",
                 date_labels = "%b %Y") +
    theme(legend.position="right") +
    geom_vline(xintercept=as.Date(controls$period_dates),linetype=2) +
    geom_label(data=plab,aes(x=pos,y=y,label=lab),size=3) 
  return(g)
}