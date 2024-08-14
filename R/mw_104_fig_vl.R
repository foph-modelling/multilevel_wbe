#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for wastewater viral load
# creation: jriou 
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_104_fig_vl = function(ww) {
  g = ww %>% 
    dplyr::filter(below_lod==0,below_loq==0) %>% 
    dplyr::group_by(ara_kt,week) %>% 
    dplyr::summarise(vl=median(vl,na.rm=TRUE),.groups="drop") %>% 
    ggplot() +
    geom_tile(aes(x=week,y=ara_kt,fill=vl)) +
    scale_y_discrete(limits=rev) +
    scale_fill_viridis_c(trans="log10",labels = function(x) format(x, scientific = TRUE)) +
    theme(axis.text = element_text(size=5)) +
    labs(x="Week",
         y="ARA",
         fill="Median viral load") +
    geom_vline(xintercept=ymd(controls$period_dates),linetype=2)
  return(g)
}