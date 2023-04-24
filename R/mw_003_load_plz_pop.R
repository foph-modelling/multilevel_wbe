#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load population at the plz level
# creation: jriou & mwagner
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_003_load_plz_pop = function() {
  
  library(pxR)
  pxR::read.px("../../02_data/geocommunes/data/")
  
  plz_pop = readxl::read_xlsx("../../02_data/geocommunes/je-f-21.03.01.xlsx",skip=5) %>% 
    dplyr::filter(!is.na(`Code commune`)) %>% 
    dplyr::select(bfs_nr=`Code commune`,
                  pop2=Habitants,
                  pop_density=`Densité de la population par km²`,
                  pop_0_19=`0-19 ans`,
                  pop_20_64=`20-64 ans`,
                  pop_65plus=`65 ans ou plus`)
  
  plz_corr = read_csv2("../../02_data/geocommunes/PLZO_CSV_WGS84/PLZO_CSV_WGS84.csv") %>% 
    dplyr::select(bfs_nr=`BFS-Nr`,
                  plz=PLZ) %>% 
    dplyr::distinct()
  
  plz_pop %>% 
    dplyr::left_join(plz_corr) %>% 
    dplyr::group_by(plz) %>% 
    dplyr::summarise
  
  
}