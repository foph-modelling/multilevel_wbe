#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for wastewater detection
# creation: jriou 
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_103_fig_detect = function(ww,detect_limit=0) {
  g = ww %>% 
    filter(!is.na(conc)) %>% 
    dplyr::mutate(detect=if_else(conc>detect_limit,"Yes","No"),.groups="drop") %>% 
    ggplot() +
    geom_tile(aes(x=date,y=ara_kt,fill=detect)) +
    scale_y_discrete(limits=rev) +
    theme(axis.text = element_text(size=5)) +
    labs(x="Day",
         y="ARA",
         fill="Detected")
  
  return(g)
}