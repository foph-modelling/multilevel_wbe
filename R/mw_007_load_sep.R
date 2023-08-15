#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: load SEP data
# creation: jriou
# init date: 2023-08-15
#:::::::::::::::::::::::::::::

mw_007_load_sep = function() {
  sep = readRDS("../../02_data/sep3/ssep3_user_geo.Rds")
}
  