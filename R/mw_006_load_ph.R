#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load public holiday data
# creation: jriou
# init date: 2023-05-15
#:::::::::::::::::::::::::::::

mw_006_load_ph = function() {
  
  years = 2020:2023
  pathURL = paste0("https://date.nager.at/api/v3/publicholidays/", years, "/CH")
  
  gatBankHol = function(X){
    
    bankHolidaysCH = fromJSON(X)
    
    bankHolidaysCH$counties[sapply(bankHolidaysCH$counties, is.null)] = "CH"
    
    bankHolidaysCH %>% select(date, counties) -> bankHolidaysCH
    
    N = nrow(bankHolidaysCH)
    max.cantons = max(sapply(bankHolidaysCH$counties, length))
    
    mat = as.data.frame(matrix(NA, ncol = max.cantons, nrow = N))
    
    for(i in 1:N){
      
      m = lengths(bankHolidaysCH$counties[i])
      mat[i, 1:m] = bankHolidaysCH$counties[[i]]
      
    }
    
    mat$date = bankHolidaysCH$date
    
    # long format
    data_long = gather(mat, name.col, canton, V1:V24, factor_key=TRUE)
    data_long$name.col = NULL
    data_long = data_long[complete.cases(data_long$canton),]
    data_long$hol = 1
    
    return(data_long)
  }
  
  holCH_list = lapply(pathURL, gatBankHol)
  holCH = do.call(rbind, holCH_list) %>% 
    as_tibble() %>% 
    dplyr::mutate(kt=substr(canton,4,5)) %>% 
    dplyr::select(-canton) %>% 
    dplyr::mutate(date=lubridate::ymd(date),
                  country=if_else(kt=="","CH",NA_character_),
                  kt=if_else(kt=="",NA_character_,kt)) 
  saveRDS(holCH, file = "data/public_holidays/holCH.rds")
}