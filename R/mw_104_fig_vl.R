#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for wastewater viral load
# creation: jriou 
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_104_fig_vl = function(ww) {

  Sys.setlocale("LC_ALL", "English")
  tt = ww %>% 
    dplyr::filter(below_lod==0,below_loq==0) %>% 
    dplyr::group_by(kt,wwtp_index,week) %>% 
    dplyr::summarise(mvl=median(vl,na.rm=TRUE),.groups="drop")
  
  ktgroup = tt %>% 
    group_by(kt) %>% 
    summarise(min=min(wwtp_index),
              max=max(wwtp_index)) %>% 
    mutate(bar=max+.5) %>% 
    filter(max<118)
  
  g =  ggplot(tt) +
    geom_tile(aes(x=week,y=wwtp_index,fill=mvl)) +
    scale_y_continuous(expand=c(0,0),breaks=c(1,seq(10,110,by=10),118)) +
    scale_x_date(expand=expansion(add=c(0,0)),
                 breaks="4 months",
                 date_labels = "%b %Y") +
    # scale_fill_viridis_c(trans="log10",labels = function(x) format(x, scientific = TRUE)) +
    scale_fill_viridis_c(trans="log",breaks=c(1e10,1e11,1e12,1e13,1e14,1e15)) +
    # theme(axis.text = element_text(size=5)) +
    labs(x="Time",
         y="WWTP",
         fill="Viral load") +
    geom_vline(xintercept=ymd(controls$period_dates),linetype=2) +
    geom_point(aes(x=as.Date("2023-12-16"),y=wwtp_index+.2,colour=kt),shape=15,size=2) +
    geom_vline(xintercept=as.Date("2023-12-08"),linewidth=.3) +
    geom_segment(data=ktgroup,aes(x=as.Date("2023-12-08"),xend=as.Date("2023-12-22"),y=bar,yend=bar),linewidth=.3) +
    # geom_point(aes(x=as.Date("2023-12-15"),y=wwtp_index),colour="black",fill="transparent",shape=22,size=2) +
    scale_colour_grey(guide="none")
  return(g)
} 
