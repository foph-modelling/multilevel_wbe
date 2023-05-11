#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: merge all and create new variables
# creation: jriou & mwagner
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_004_merge_all = function(ww,ms,pd) {
  # merge all
  out = ww %>% 
    dplyr::left_join(ms, by = join_by(ara_id, date)) %>% 
    dplyr::left_join(pd, by="ara_id")
  
  # create vl with new population
  out = out %>% 
    dplyr::mutate(vl=(conc*flow*1e3)/(pop_total/1e5)) 

  return(out)
}