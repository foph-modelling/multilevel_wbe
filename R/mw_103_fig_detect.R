#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for wastewater detection
# creation: jriou 
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_103_fig_detect = function(ww) {
  g = ww %>% 
    dplyr::mutate(detect=case_when(.default = "Detected",
                                   below_loq==1 ~ "Under LOQ",
                                   below_lod==1 ~ "Under LOD"),.groups="drop") %>% 
    ggplot() +
    geom_tile(aes(x=date,y=ara_kt,fill=detect)) +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values=c("grey","steelblue","firebrick")) +
    labs(x="Date",y="ARA",fill="Wastewater SARS-CoV-2")

  return(g)
}