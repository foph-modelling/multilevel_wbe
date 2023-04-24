#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for wastewater concentration
# creation: jriou 
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_102_fig_conc = function(ww,lower_limit=0) {
  g = ww %>% 
    dplyr::group_by(ara_kt,week) %>% 
    dplyr::summarise(conc=mean(conc,na.rm=TRUE),.groups="drop") %>% 
    dplyr::filter(!is.na(conc),conc>lower_limit) %>%
    ggplot() +
    geom_tile(aes(x=week,y=ara_kt,fill=conc)) +
    scale_y_discrete(limits=rev) +
    scale_fill_viridis_c(trans="log10") +
    theme(axis.text = element_text(size=5)) +
    labs(x="Week",
         y="ARA",
         fill="Mean log concentration")
  
  return(g)
}