#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load wastewater data
# creation: jriou
# init date: 2023-03-08
#:::::::::::::::::::::::::::::

mw_001_load_ww = function() {
  
  # setup paths ----
  ldate = fs::path(paste0(controls$data_path,controls$data_date))
<<<<<<< HEAD
  
  # load raw data ----
  raw_ww = readRDS(fs::path(paste0(ldate,"_ww_data/COVID19Wastewater_vl.rds")))
=======
  data_context_url = "https://www.covid19.admin.ch/api/data/context"
  data_context_json = httr::GET(data_context_url) %>% 
    jsonlite::parse_json()
  ww_url = data_context_json$sources$individual$csv$wasteWater$viralLoad
  
  # load raw data ----
  raw_ww = data.table::fread(ww_url) %>% 
    dplyr::as_tibble()
  # raw_ww = readRDS(fs::path(paste0(ldate,"_ww_data/COVID19Wastewater_vl.rds")))
>>>>>>> gh-pages
  ara_supp = readRDS(fs::path(paste0(ldate,"_ww_data/ara_supp_tbl.rds"))) %>% 
    dplyr::mutate(ara_n=row_number(),
                  lab_n=as.factor(as.numeric(as.factor(lab))))
  ara_methods = readRDS(fs::path(paste0(ldate,"_ww_data/ara_method_change_dates.rds")))
<<<<<<< HEAD
=======
  # periods
  cut_dates = lubridate::ymd(c(min(raw_ww$date)-1,controls$period_dates),max(raw_ww$date)+1)
  # NUTS
  nuts_regions = data.frame(kt = c("VD","VS","GE","BE","FR","SO","NE","JU","BS","BL","AG","ZH","GL",
                               "SH","AR","AI","SG","GR","TG","LU","UR","SZ","OW","NW","ZG","TI"),
                    NUTS2 = c(rep(c(1:7),c(3,5,3,1,7,6,1))),
                    NUTS2_name = rep(c("Lake Geneva","Mittelland","Northwest","Zurich","Eastern","Central","Ticino"),c(3,5,3,1,7,6,1)))
>>>>>>> gh-pages
  
  # initial management ----
  ww = raw_ww %>% 
    dplyr::select(ara_id=geoRegion,
                  ara_name=name,
                  kt=canton,
                  date,
                  version,
<<<<<<< HEAD
                  flow,conc,vl) %>% 
    # geographic information
    dplyr::mutate(ara_kt=paste0(kt," / ",ara_name)) %>% 
    # time information
    dplyr::mutate(isoweek=ISOweek::date2ISOweek(date),
                  week=ISOweek::ISOweek2date(paste0(substr(isoweek,1,9),"1"))) %>% 
    dplyr::select(-isoweek) %>% 
=======
                  flow,conc,vl,
                  below_loq,loq,
                  method_date_from,method_date_to) %>% 
    # geographic information
    dplyr::mutate(ara_kt=paste0(kt," / ",ara_name)) %>% 
    # time information
    dplyr::mutate(date=lubridate::ymd(date),
                  isoweek=ISOweek::date2ISOweek(date), # get week in ISO format
                  week=ISOweek::ISOweek2date(paste0(substr(isoweek,1,9),"4")), # attribute week to the date of the Thursday (mid-week)
                  ara_id=as.character(ara_id)) %>% 
>>>>>>> gh-pages
    # add information about labs
    dplyr::left_join(ara_supp, by = join_by(ara_id, ara_name, kt)) %>% 
    dplyr::arrange(ara_n,ara_name,kt,date) %>% 
    dplyr::relocate(ara_n) %>% 
    # add indicator for methods change
    dplyr::left_join(ara_methods, by = join_by(ara_id, ara_name)) %>% 
    dplyr::group_by(ara_id,kt,lab) %>% 
    dplyr::mutate(method=if_else(date<method_date_change,0,1),
<<<<<<< HEAD
                  method=if_else(is.na(method),0,method)) %>% 
    # add indicator for missing/
    dplyr::mutate(measurement=if_else(is.na(conc),0,1)) %>%
    # remove Liechtenstein (strange values and not CH) and Ramsen (mostly in Germany)
    dplyr::filter(ara_id!="100000", ara_id!="296300") %>%
    dplyr::ungroup()
=======
                  method=if_else(is.na(method),0,method)) %>%  # within an ARA, 1 indicates a new method compared to baseline (only 1 method change observed)
    # deal with LOD and LOQ 
    # TODO: flag values below LOQ in other ARAs (using the lowest LOQ?)
    dplyr::mutate(
      measurement=if_else(is.na(conc),0,1),
                  below_lod=conc==0) %>%
    # remove Liechtenstein (not CH), Ramsen (mostly in Germany) and Burgdorf (issues with initial method, reintegrated after method change on 2022-08-23)
    dplyr::filter(ara_id!="100000", 
                  ara_id!="296300",
                  !(ara_id==40100 & date<"2022-08-23")) %>%
    # create periods from controls 
    dplyr::mutate(period=cut(date,breaks=cut_dates),
                  period=as.factor(as.numeric(as.factor(period)))) %>% 
    # identify weekends and holidays
    # TODO: put holidays as well
    mutate(day=as.numeric(date)-min(as.numeric(date)),
           weekend=if_else(lubridate::wday(date,week_start=1)>=5,1,0)) %>% 
    # create regions
    dplyr::left_join(nuts_regions,by="kt")
    
  
>>>>>>> gh-pages
  return(ww)
}