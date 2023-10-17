#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: entry point
# creation: jriou
# init date: 2023-03-08
#:::::::::::::::::::::::::::::


# Block 0: controls and set-up --------------------------------------------

analysis_date = "2023-09-29"
data_date = "2023-09-29" # name of the repertory with fixed data (up-to-date wastewater data is downloaded directly at every run)
data_path = "../../02_data/wastewater/"
period_dates = c("2022-05-16","2022-09-05","2023-01-02","2023-07-03") # set at the lowest points between waves, on Mondays so weeks are not cut
controls = list(update_data=TRUE, # set to TRUE before sourcing to update the data
                rerun_models=TRUE, # only applies to large models
                compute_cv=FALSE, 
                data_date=data_date,
                analysis_date=analysis_date,
                data_path=data_path,
                period_dates=period_dates,
                savepoint=paste0("savepoints/savepoint_",analysis_date,"/")) # create new repertory whenever analysis_date changes
if(FALSE) saveRDS(controls, paste0(controls$savepoint,"controls.rds"))
source("R/setup.R")

# Block 1: data prep ------------------------------------------------------

if(controls$update_data) {
  # load and prep wastewater data
  ww0 = mw_001_load_ww()
  # load and prep corresponding reported tests, cases and hospits
  ms0 = mw_002_load_ms()
  # load and prep population data by ARA # TODO: add population by hectare
  pd0 = mw_003_load_plz_pop()
  # load SEP
  se0 = mw_007_load_sep()
  # merge and last prep
  ww1 = mw_004_merge_all(ww0,ms0,pd0)
  # load shape files
  shapes = mw_005_load_shp()
  # save
  saveRDS(ww1,file=fs::path(controls$savepoint,"ww1.rds"))
  saveRDS(shapes,file=fs::path(controls$savepoint,"shapes.rds"))
}


# Block 2: analysis ----------------------------------------------------

# load data for development (to be removed)
if(FALSE) {
  ww1 = readRDS(fs::path(controls$savepoint,"ww1.rds"))
  shapes = readRDS(fs::path(controls$savepoint,"shapes.rds"))
}

# data description
rmarkdown::render("R/mw_200_data_description.R",
                  params=list(controls=controls),
                  output_file=file.path("../",controls$savepoint,"data_description.html"),clean=FALSE)

# model development
rmarkdown::render("R/mw_201_model_dev_A1.R",
                  params=list(controls=controls),
                  output_file=file.path("../",controls$savepoint,"model_dev_A1.html"),clean=FALSE)
rmarkdown::render("R/mw_201_model_dev_A2.R",
                  params=list(controls=controls),
                  output_file=file.path("../",controls$savepoint,"model_dev_A2.html"),clean=FALSE)
rmarkdown::render("R/mw_201_model_dev_A3.R",
                  params=list(controls=controls),
                  output_file=file.path("../",controls$savepoint,"model_dev_A3.html"),clean=FALSE)

