#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for missing and available wastewater data
# creation: jriou 
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_101_fig_missing = function(ww,identify_eawag=FALSE) {
  g = ww %>% 
    ggplot() +
    geom_point(aes(x=date,y=ara_kt,colour=kt),size=.5) +
    scale_y_discrete(limits=rev) +
    scale_colour_discrete(guide="none") +
    theme(axis.text = element_text(size=5)) +
    labs(x="Date",y="ARA")
  if(identify_eawag) {
    g = ww %>% 
      mutate(eawag=ifelse(lab=="EAWAG","EAWAG","other")) %>% 
      ggplot() +
      geom_point(aes(x=date,y=ara_kt,fill=kt,colour=eawag),shape=21,size=.8,stroke=.1) +
      scale_y_discrete(limits=rev) +
      scale_colour_manual(values=c("black","transparent")) +
      scale_fill_discrete(guide="none") +
      theme(axis.text = element_text(size=5)) +
      labs(x="Date",y="ARA",colour="Laboratory")
  }
  
  return(g)
}