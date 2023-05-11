#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load mandatory surveillance data
# creation: jriou & mwagner
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_002_load_ms = function() {
  
  # setup paths ----
  ldate = paste0(controls$data_path,controls$data_date)
  
  # load raw data ----
  plz_ara = readRDS(paste0(ldate,"_ww_data/plz_ara.rds")) %>% 
    dplyr::select(ara_id,plz,per_plz)
  ncov2019_falldetail_cases_plz = readRDS(paste0(ldate,"_ww_data/ncov2019_falldetail_cases_plz.rds"))
  ncov2019_labortests_plz = readRDS(paste0(ldate,"_ww_data/ncov2019_labortests_plz.rds"))
  ncov2019_neg_labortests_plz = readRDS(paste0(ldate,"_ww_data/ncov2019_neg_labortests_plz.rds"))
  
  # aggregating function ----
  aggr_ara <- function(data, plz_ara, ...) {
    data_count <-
      data %>%
      dplyr::left_join(plz_ara, by = c("plz_pat" = "plz"),relationship = "many-to-many") %>%
      dplyr::filter(!is.na(ara_id)) %>%
      dplyr::group_by(ara_id, ...) %>%
      dplyr::summarise(count_per_ara = sum(per_plz/100),.groups="drop") %>%
      dplyr::ungroup()
    return(data_count)
  }
  
  # cases by ara by date ----
  ms_reported_cases = ncov2019_falldetail_cases_plz %>%
    aggr_ara(plz_ara, fall_dt) %>% 
    dplyr::select(ara_id,date=fall_dt,reported_cases=count_per_ara)
  
  # hospits by ara by date ----
  ms_reported_hospit = ncov2019_falldetail_cases_plz %>%
    aggr_ara(plz_ara, pt_hospdatin) %>%
    filter(!is.na(pt_hospdatin)) %>% 
    dplyr::select(ara_id,date=pt_hospdatin,reported_hospit=count_per_ara)
  
  # tests by ara by date ----
  ms_reported_tests = bind_rows(ncov2019_labortests_plz, ncov2019_neg_labortests_plz) %>%
    aggr_ara(plz_ara, test_dt) %>% 
    dplyr::select(ara_id,date=test_dt,reported_tests=count_per_ara)
  
  # merge
  ms_reported = tidyr::expand_grid(ara_id=unique(ms_reported_cases$ara_id),
                                   date=seq.Date(from=min(ms_reported_cases$date),
                                                 to=max(ms_reported_cases$date),
                                                 by=1)) %>% 
    dplyr::left_join(ms_reported_tests,by = join_by(ara_id, date)) %>% 
    dplyr::left_join(ms_reported_cases,by = join_by(ara_id, date)) %>% 
    dplyr::left_join(ms_reported_hospit,by = join_by(ara_id, date)) %>% 
    tidyr::replace_na(list(reported_cases=0,reported_hospit=0,reported_tests=0))
    
  return(ms_reported)
}