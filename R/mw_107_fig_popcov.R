#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for population covered
# creation: jriou 
# init date: 2023-11-04
#:::::::::::::::::::::::::::::

mw_107_fig_popcov = function(ww) {
  # gather info
  tt = ww %>% 
    dplyr::group_by(week,ara_name) %>% 
    dplyr::summarise(pop=max(pop),.groups="drop_last") %>% 
    dplyr::summarise(pop=sum(pop),
                     pop_cov=sum(pop)/8.7e6)
  
  # plot fig
  g = ggplot(tt) +
    geom_line(aes(x=week,y=pop_cov),colour=cust_cols[3]) +
    scale_y_continuous(label=scales::percent,limits=c(0,1)) +
    labs(x="Date",y="Population covered")
  
  return(g)
}