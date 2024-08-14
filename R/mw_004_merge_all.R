#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: merge all and create new variables
# creation: jriou & mwagner
# init date: 2023-04-17
#:::::::::::::::::::::::::::::


mw_004_merge_all = function(ww,ms,pd,pc, complete_dates=FALSE) {
  
  if(complete_dates==TRUE){
    ww = ww %>% complete(ara_id, day)
  }

  # merge all
  se = se %>% dplyr::select(-ara_n)
  out = ww %>% 
    dplyr::left_join(ms, by = join_by(ara_id, date)) %>% 
    dplyr::left_join(pd, by="ara_id") %>% 
    dplyr::left_join(pc, by="ara_id") %>% 
    dplyr::mutate(ara_n=as.character(ara_n)) %>% 
    dplyr::left_join(se, by="ara_id") %>% 
    dplyr::left_join(co, by="ara_id") %>% 

    dplyr::ungroup()
  
  # create vl with new population
  out = out %>% 
    dplyr::mutate(vl=(conc*flow*1e3)/(pop_total/1e5)) 
  
  # create employment_per_population
  # out = out %>% 
    # dplyr::mutate(employment_factor=catch_emp_pop/pop_total)

  if(FALSE) {
    out %>% 
      group_by(ara_name) %>% 
      summarise(pop=max(pop),pop_total=max(pop_total)) %>% 
      mutate(ratio=pop_total/pop,
             absratio=abs(1-ratio)) %>% 
      arrange(-absratio) %>% 
      view()
    
    out %>% 
      ungroup() %>% 
      summarise(pop=sum(pop),
                pop2=sum(pop_total,na.rm=TRUE))
  }
  return(out)
}
