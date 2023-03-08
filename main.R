#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: entry point
# creation: jriou
# init date: 2023-03-08
#:::::::::::::::::::::::::::::

# Block 0: set-up ----
analysis_date = "2023-03-08"
controls = list(update_data=TRUE,
                savepoint=paste0("savepoint_",analysis_date,"/"))
source("R/setup.R")

# Block 1: load data ----
if(controls$update_data) {
  mw0 = mw_001_load_data()
}

# Block 2: data management ----

# Block 3: description ----

# Block 4: models ----

