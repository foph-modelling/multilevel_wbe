#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load wastewater data
# creation: jriou
# init date: 2023-03-08
#:::::::::::::::::::::::::::::

mw_001_load_ww = function() {
  
  # setup paths ----
  ldate = fs::path(paste0(controls$data_path,controls$data_date))
  
  # load raw data ----
  raw_ww = readRDS(fs::path(paste0(ldate,"_ww_data/COVID19Wastewater_vl.rds")))
  ara_supp = readRDS(fs::path(paste0(ldate,"_ww_data/ara_supp_tbl.rds"))) %>% 
    dplyr::mutate(ara_n=row_number(),
                  lab_n=as.factor(as.numeric(as.factor(lab))))
  ara_methods = readRDS(fs::path(paste0(ldate,"_ww_data/ara_method_change_dates.rds")))
  
  # initial management ----
  ww = raw_ww %>% 
    dplyr::select(ara_id=geoRegion,
                  ara_name=name,
                  kt=canton,
                  date,
                  version,
                  flow,conc,vl) %>% 
    # geographic information
    dplyr::mutate(ara_kt=paste0(kt," / ",ara_name)) %>% 
    # time information
    dplyr::mutate(isoweek=ISOweek::date2ISOweek(date),
                  week=ISOweek::ISOweek2date(paste0(substr(isoweek,1,9),"1"))) %>% 
    dplyr::select(-isoweek) %>% 
    # add information about labs
    dplyr::left_join(ara_supp, by = join_by(ara_id, ara_name, kt)) %>% 
    dplyr::arrange(ara_n,ara_name,kt,date) %>% 
    dplyr::relocate(ara_n) %>% 
    # add indicator for methods change
    dplyr::left_join(ara_methods, by = join_by(ara_id, ara_name)) %>% 
    dplyr::group_by(ara_id,kt,lab) %>% 
    dplyr::mutate(method=if_else(date<method_date_change,0,1),
                  method=if_else(is.na(method),0,method)) %>% 
    # add indicator for missing/
    dplyr::mutate(measurement=if_else(is.na(conc),0,1)) %>%
    # remove Liechtenstein (strange values and not CH) and Ramsen (mostly in Germany)
    dplyr::filter(ara_id!="100000", ara_id!="296300") %>%
    dplyr::ungroup()
  return(ww)
}