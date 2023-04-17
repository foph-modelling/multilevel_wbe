#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: entry point
# creation: jriou
# init date: 2023-03-08
#:::::::::::::::::::::::::::::

# Block 0: set-up ----
data_date = "2023-03-14"
analysis_date = "2023-04-17"
data_path = "../../02_data/wastewater/"
controls = list(update_data=TRUE,
                data_date=data_date,
                data_path=data_path,
                savepoint=paste0("savepoint_",analysis_date,"/"))
source("R/setup.R")

# Block 1: data prep ----
if(controls$update_data) {
  # load and prep wastewater data
  ww0 = mw_001_load_ww()
  # load, prep and merge corresponding reported tests, cases and hospits
  ms0 = mw_002_load_ms()
  ww1 = left_join(ww0,ms0,by = join_by(ara_id, date))
}

# Block 2: description ----
display_items = list(
  
)

# Block 3: models ----

