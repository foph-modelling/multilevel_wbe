#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for wastewater detection
# creation: jriou 
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_103_fig_detect = function(ww,detect_limit=0) {
  g = ww %>% 
    filter(!is.na(conc)) %>% 
    dplyr::mutate(detect=case_when(.default = "Detected",
                                   below_loq==TRUE ~ "Under LOQ",
                                   below_lod==TRUE ~ "Under LOD"),.groups="drop") %>% 
    ggplot() +
    geom_tile(aes(x=date,y=ara_kt,fill=detect)) +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values=c("grey","steelblue","firebrick"))

  return(g)
}