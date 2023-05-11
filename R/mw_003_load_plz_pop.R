#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load population at the plz level
# creation: jriou & mwagner
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_003_load_plz_pop = function() {
  
<<<<<<< HEAD
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
  
  
=======
  # setup paths ----
  ldate = paste0(controls$data_path,controls$data_date)
  
  # load raw data from BFS https://www.bfs.admin.ch/asset/de/su-d-01.02.03.07 ----
  plz_pop = readxl::read_xlsx("data/su-d-01.02.03.07.xlsx",skip=2) 
  
  # load correspondence between PLZ and ARA
  plz_ara = readRDS(paste0(ldate,"_ww_data/plz_ara.rds")) %>% 
    dplyr::select(ara_id,plz,per_plz)
  
  # reformat ----
  plz_pop = plz_pop %>% 
    dplyr::rename(plz=1) %>% 
    dplyr::filter(plz!="Schweiz",plz!="100000",!is.na(Total)) %>% 
    dplyr::mutate(plz=as.numeric(plz),
                  pop_0_19=`0-4`+`5-9`+`10-14`+`15-19`,
                  pop_20_64=`20-24`+`25-29`+`30-34`+`35-39`+`40-44`+`45-49`+`50-54`+`55-59`+`60-64`,
                  pop_65plus=`65-69`+`70-74`+`75-79`+`80-84`+`85-89`+`90 und mehr`) %>% 
    # dplyr::filter(plz!="Schweiz",plz!="100000",!is.na(plz)) %>% 
    dplyr::select(plz,pop_total=Total,pop_0_19,pop_20_64,pop_65plus)
  
  # aggregate from PLZ to ARA
  ara_pop = plz_ara %>% 
    dplyr::left_join(plz_pop,by="plz") %>% 
    dplyr::group_by(ara_id) %>% 
    dplyr::summarise(pop_total=sum(pop_total*per_plz/100),
                     pop_0_19=sum(pop_0_19*per_plz/100),
                     pop_20_64=sum(pop_20_64*per_plz/100),
                     pop_65plus=sum(pop_65plus*per_plz/100))%>% 
    dplyr::filter(ara_id !="100000") 
  
  return(ara_pop)
>>>>>>> gh-pages
}