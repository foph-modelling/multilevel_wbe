#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: initialize R session
# creation: jriou
# init date: 2023-03-08
#:::::::::::::::::::::::::::::


# load libraries ----
pacman::p_load(Hmisc,
               data.table,
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
               scales, 
               units)

# set language ----

Sys.setenv(LANG = "en")

# source functions ----
fili = dir(path = "R",pattern="mw_[01456789]",full.names = TRUE)
lapply(X = fili, FUN = function(x) {source(paste0(x), echo=FALSE)})

# create savepoint repertory if not existing ----
dir.create(file.path(controls$savepoint), showWarnings = FALSE)
# saveRDS(controls,file=fs::path(controls$savepoint,"controls.rds"))

# common ----

# small custom functions ----
qsum_range = function(x) {
  if(median(x,na.rm=TRUE)>1e6) {
  r = paste0(formatC(median(x,na.rm=TRUE), format="e",big.mark=",", digits=1),
             " (range: ",
             formatC(min(x,na.rm=TRUE), format="e", big.mark=",", digits=1),
             " to ",
             formatC(max(x,na.rm=TRUE), format="e", big.mark=",", digits=1),
             ")") 
  } else if (median(x,na.rm=TRUE)>1) {
    r = paste0(formatC(median(x,na.rm=TRUE), format="f",big.mark=",", digits=0),
               " (range: ",
               formatC(min(x,na.rm=TRUE), format="f", big.mark=",", digits=0),
               " to ",
               formatC(max(x,na.rm=TRUE), format="f", big.mark=",", digits=0),
               ")") 
  } else {
    r = paste0(formatC(median(x,na.rm=TRUE), format="f",big.mark=",", digits=2),
               " (range: ",
               formatC(min(x,na.rm=TRUE), format="f", big.mark=",", digits=2),
               " to ",
               formatC(max(x,na.rm=TRUE), format="f", big.mark=",", digits=2),
               ")") 
  }
}
# aesthetics ----
theme_set(theme_bw())
cust_cols = c("dodgerblue","firebrick2","goldenrod2")