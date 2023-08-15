#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: initialize R session
# creation: jriou
# init date: 2023-03-08
#:::::::::::::::::::::::::::::


# load libraries ----
pacman::p_load(data.table,
               tidyverse,
               lubridate,
               ISOweek,
               INLA,
               inlabru,
               sf,
               splines,
               cowplot,
               #flextable,
               spdep,
               jsonlite,
               scales)

# source functions ----
fili = dir(pattern="mw_[013456789]")
lapply(X = fili, FUN = function(x) {source(paste0(x), echo=FALSE)})

# create savepoint repertory if not existing ----
dir.create(file.path("./", controls$savepoint), showWarnings = FALSE)
# saveRDS(controls,file=fs::path(controls$savepoint,"controls.rds"))

# common ----

# small custom functions ----
qsum_range = function(x) {
    r = paste0(formatC(median(x,na.rm=TRUE), format="g",big.mark=",", digits=0),
               " (range: ",
               formatC(min(x,na.rm=TRUE), format="g", big.mark=",", digits=0),
               " to ",
               formatC(max(x,na.rm=TRUE), format="g", big.mark=",", digits=0),
               ")")
}

# aesthetics ----
theme_set(theme_bw())
cust_cols = c("dodgerblue","firebrick2","goldenrod2")
