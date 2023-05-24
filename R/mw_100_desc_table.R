#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create descriptive table for wastewater data
# creation: jriou 
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_100_desc_table = function(ww, ...) {
  tab = ww %>% 
    dplyr::filter(!is.na(vl)) %>% 
    dplyr::group_by(...) %>% 
    dplyr::summarise(
      `Number of ARAs`=length(unique(ara_kt)),
      `Number of laboratories`=length(unique(lab)),
      `Number of laboratory methods`=length(unique(paste0(lab,method))),
      `Number of measurements`=sum(!is.na(vl)),
      `Measurements below LOQ`=sum(below_loq==1),
      `Measurements below LOD`=sum(below_lod==1),
      `First`=min(date,na.rm=TRUE),
      `Last`=max(date,na.rm=TRUE),
      `Median viral concentration [gc/L]`=qsum_range(conc),
      `Median flow [m3/day]`=qsum_range(flow),
      `Median viral load [gc/day/100,000]`=qsum_range(vl))
  return(tab)
}
