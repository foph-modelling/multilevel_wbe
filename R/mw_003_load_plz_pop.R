#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load population at the plz level
# creation: jriou & mwagner
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_003_load_plz_pop = function() {
  
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
    dplyr::select(plz,pop_total=Total,pop_0_19,pop_20_64,pop_65plus)
  
  # aggregate from PLZ to ARA
  ara_pop = plz_ara %>% 
    dplyr::left_join(plz_pop,by="plz") %>% 
    dplyr::group_by(ara_id) %>%
    dplyr::filter(!is.na(pop_total)) %>% 
    dplyr::summarise(pop_total=sum(pop_total*per_plz/100),
                     pop_0_19=sum(pop_0_19*per_plz/100),
                     pop_20_64=sum(pop_20_64*per_plz/100),
                     pop_65plus=sum(pop_65plus*per_plz/100))%>% 
    dplyr::filter(ara_id !="100000") 
  
  return(ara_pop)
}