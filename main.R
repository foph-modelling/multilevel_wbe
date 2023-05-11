#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: entry point
# creation: jriou
# init date: 2023-03-08
#:::::::::::::::::::::::::::::


<<<<<<< HEAD
# Block 0: set-up ---------------------------------------------------------

data_date = "2023-03-14"
analysis_date = "2023-04-17"
data_path = "../../02_data/wastewater/"
controls = list(update_data=FALSE,
                data_date=data_date,
                analysis_date=analysis_date,
                data_path=data_path,
                savepoint=paste0("savepoints/savepoint_",analysis_date,"/"))
=======
# Block 0: controls and set-up --------------------------------------------

analysis_date = "2023-05-09"
data_date = "2023-03-14" # name of the repertory with fixed data (up-to-date wastewater data is downloaded directly at every run)
data_path = "../../02_data/wastewater/"
period_dates = c("2022-05-16","2022-09-05","2023-01-02") # set at the lowest points between waves, on Mondays so weeks are not cut
controls = list(update_data=FALSE, # set to TRUE before sourcing to update the data
                data_date=data_date,
                analysis_date=analysis_date,
                data_path=data_path,
                period_dates=period_dates,
                savepoint=paste0("savepoints/savepoint_",analysis_date,"/")) # create new repertory whenever analysis_date changes
>>>>>>> gh-pages
source("R/setup.R")

# Block 1: data prep ------------------------------------------------------

if(controls$update_data) {
  # load and prep wastewater data
  ww0 = mw_001_load_ww()
<<<<<<< HEAD
  # load, prep and merge corresponding reported tests, cases and hospits
  ms0 = mw_002_load_ms()
  ww1 = left_join(ww0,ms0,by = join_by(ara_id, date))
  saveRDS(ww1,file=fs::path(controls$savepoint,"ww1.rds"))
=======
  # load and prep corresponding reported tests, cases and hospits
  ms0 = mw_002_load_ms()
  # load and prep population data by ARA
  pd0 = mw_003_load_plz_pop()
  # merge and last prep
  ww1 = mw_004_merge_all(ww0,ms0,pd0)
  # load shape files
  shapes = mw_005_load_shp()
  # save
  saveRDS(ww1,file=fs::path(controls$savepoint,"ww1.rds"))
  saveRDS(shapes,file=fs::path(controls$savepoint,"shapes.rds"))
>>>>>>> gh-pages
}


# Block 2: analysis ----------------------------------------------------

<<<<<<< HEAD
=======
if(FALSE) {
  ww1 = readRDS(fs::path(controls$savepoint,"ww1.rds"))
  shapes = readRDS(fs::path(controls$savepoint,"shapes.rds"))
}

>>>>>>> gh-pages
# data description
rmarkdown::render("R/mw_200_data_description.R",
                  params=list(controls=controls),
                  output_file="../reports/data_description.html")

# model development
<<<<<<< HEAD
rmarkdown::render("mw_201_model_dev_A.R",
=======
rmarkdown::render("R/mw_201_model_dev_A.R",
>>>>>>> gh-pages
                  params=list(controls=controls),
                  output_file="../reports/model_dev_A.html")



