#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for wastewater viral load
# creation: jriou 
# init date: 2023-08-07
#:::::::::::::::::::::::::::::

mw_106_fig_vl_time = function(ww) {
  g = ww %>% 
    dplyr::filter(below_lod==0,below_loq==0) %>% 
    # dplyr::group_by(ara_kt,week) %>% 
    # dplyr::summarise(vl=mean(vl,na.rm=TRUE),.groups="drop") %>% 
    ggplot() +
    geom_line(aes(x=date,y=vl,colour=vl,group=ara_id),alpha=.5) +
    scale_colour_viridis_c(trans="log",guide="none") +
    scale_y_continuous(trans="log",
                       breaks = trans_breaks("log", function(x) exp(x)),
                       labels = trans_format("log", math_format(e^.x))) +  
    labs(x="Week",
         y="Viral load",
         fill="Mean viral load")
  return(g)
}