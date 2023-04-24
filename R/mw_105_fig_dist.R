#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for wastewater viral load distribution
# creation: jriou 
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_105_fig_dist = function(ww,lower_limit=0) {
  g = ww %>% 
    dplyr::filter(!is.na(vl),vl>lower_limit) %>%
    ggplot() +
    geom_histogram(aes(x=log10(conc)))
    
  return(g)
}