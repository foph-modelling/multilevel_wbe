#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load wastewater data
# creation: jriou
# init date: 2023-03-08
#:::::::::::::::::::::::::::::

mw_001_load_ww = function() {
  
  # setup paths ----
  ldate = fs::path(paste0(controls$data_path,controls$data_date))
  data_context_url = "https://www.covid19.admin.ch/api/data/context"
  data_context_json = httr::GET(data_context_url) %>% 
    jsonlite::parse_json()
  ww_url = data_context_json$sources$individual$csv$wasteWater$viralLoad
  
  # load raw data ----
  raw_ww = data.table::fread(ww_url) %>% 
    dplyr::as_tibble()
  # raw_ww = readRDS(fs::path(paste0(ldate,"_ww_data/COVID19Wastewater_vl.rds")))
  ara_supp = readRDS(fs::path(paste0(ldate,"_ww_data/ara_supp_tbl.rds"))) %>% 
    dplyr::mutate(ara_n=row_number(),
                  lab_n=as.factor(as.numeric(as.factor(lab))))
  ara_methods = readRDS(fs::path(paste0(ldate,"_ww_data/ara_method_change_dates.rds")))
  # periods
  cut_dates = lubridate::ymd(c(min(raw_ww$date)-1,controls$period_dates),max(raw_ww$date)+1)
  # NUTS
  nuts_regions = data.frame(kt = c("VD","VS","GE","BE","FR","SO","NE","JU","BS","BL","AG","ZH","GL",
                                   "SH","AR","AI","SG","GR","TG","LU","UR","SZ","OW","NW","ZG","TI"),
                            NUTS2 = c(rep(c(1:7),c(3,5,3,1,7,6,1))),
                            NUTS2_name = rep(c("Lake Geneva","Mittelland","Northwest","Zurich","Eastern","Central","Ticino"),c(3,5,3,1,7,6,1)))
  # national and cantonal holidays
  mw_006_load_ph() # download updated list of national and cantonal holidays from date.nager.at
  hol = readRDS(file = "data/public_holidays/holCH.rds")
  hol_CH = hol %>% 
    dplyr::filter(country=="CH") %>% 
    dplyr::select(date,hol1=hol)
  hol_kt = hol %>% 
    dplyr::filter(is.na(country)) %>% 
    dplyr::select(kt,date,hol2=hol)
  
  # initial management ----
  ww = raw_ww %>% 
    dplyr::select(ara_id=geoRegion,
                  ara_name=name,
                  kt=canton,
                  date,
                  version,
                  flow,conc,vl,
                  below_loq,loq,
                  method_date_from,method_date_to) %>% 
    # remove missing measurement
    dplyr::filter(!is.na(conc)) %>% 
    # geographic information
    dplyr::mutate(ara_kt=paste0(kt," / ",ara_name)) %>% 
    # time information
    dplyr::mutate(date=lubridate::ymd(date),
                  isoweek=ISOweek::date2ISOweek(date), # get week in ISO format
                  week=ISOweek::ISOweek2date(paste0(substr(isoweek,1,9),"4")), # attribute week to the date of the Thursday (mid-week)
                  ara_id=as.character(ara_id)) %>% 
    dplyr::filter(date<=lubridate::ymd(controls$analysis_date)) %>%
    # add information about labs
    dplyr::left_join(ara_supp, by = join_by(ara_id, ara_name, kt)) %>% 
    dplyr::arrange(ara_n,ara_name,kt,date) %>% 
    dplyr::relocate(ara_n) %>% 
    # add indicator for methods change
    dplyr::left_join(ara_methods, by = join_by(ara_id, ara_name)) %>% 
    dplyr::group_by(ara_id,kt,lab) %>% 
    dplyr::mutate(method=if_else(date<method_date_change,0,1),
                  method=if_else(is.na(method),0,method)) %>%  # within an ARA, 1 indicates a new method compared to baseline (only 1 method change observed)
    dplyr::ungroup() %>% 
    # deal with LOD and LOQ 
    dplyr::mutate(loq=if_else(loq==27.4,NA_real_,loq),  # remove mistake with LOQ at 27.4 in Aire
                  loq=if_else(is.na(loq),700,loq), # use a LOQ of 700 if missing (as in EAWAG)
                  below_loq=below_loq | (conc<=loq & conc>0),
                  conc=if_else(below_loq,loq,conc),
                  below_loq=as.factor(if_else(below_loq==TRUE,1,0)),
                  below_lod=as.factor(if_else(conc==0,1,0,0))) %>%
    dplyr::relocate(below_lod, .after=vl) %>% 
    # remove strange values in canton NE around April-May 2022
    dplyr::filter(!(kt=="NE" & conc<1)) %>% 
    # remove Liechtenstein (not CH), Ramsen (mostly in Germany) and Burgdorf (issues with initial method, reintegrated after method change on 2022-08-23)
    dplyr::filter(ara_id!="100000", 
                  ara_id!="296300",
                  !(ara_id==40100 & date<"2022-08-23")) %>%
    # create periods from controls 
    dplyr::mutate(period=cut(date,breaks=cut_dates),
                  period=as.factor(as.numeric(as.factor(period)))) %>% 
    # identify weekends, national and cantonal holidays
    dplyr::mutate(day=as.numeric(date)-min(as.numeric(date)),
                  weekend=if_else(lubridate::wday(date,week_start=1)>=5,1,0)) %>% 
    dplyr::left_join(hol_CH,by="date") %>% 
    dplyr::left_join(hol_kt,by=c("kt","date")) %>% 
    dplyr::mutate(hol=if_else(hol1==1 | hol2==1,1,0,0)) %>% 
    dplyr::select(-hol1,-hol2) %>% 
    # create regions
    dplyr::left_join(nuts_regions,by="kt")

  return(ww)
}