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
      `Total measurements`=sum(!is.na(vl)),
      `Measurements below LOQ`=sum(below_loq==1),
      `Measurements below LOD`=sum(below_lod==1),
      `Date of first measurement`=min(date,na.rm=TRUE),
      `Date of last measurement`=max(date,na.rm=TRUE),
      `Viral concentration [gc/L]`=qsum_range(conc),
      `Wastewater flow [m3/day]`=qsum_range(flow),
      `Viral load [gc/day/100,000]`=qsum_range(vl),
      `Population covered`=qsum_range(pop_total),
      # `Population density`=qsum_range(pop_dens),
      `Proportion of population under 20`=qsum_range(prop_under_20),
      `Proportion of population over 65`=qsum_range(prop_over_65),
      `Proportion of population non-Swiss`=qsum_range(prop_non_ch_eu),
      `Median of Swiss index of socio-economic position (SEP)`=qsum_range(ssep3_med),
      # `Standard deviation of Swiss index of socio-economic position (SEP)`=qsum_range(ssep3_sd),
      `Median employment factor`=qsum_range(employment_factor)
      )
  return(tab)
}
