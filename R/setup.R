#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: initialize R session
# creation: jriou
# init date: 2023-03-08
#:::::::::::::::::::::::::::::


# load libraries ----
library(tidyverse)
library(lubridate)
library(ISOweek)
library(INLA)
library(sf)

# set paths ----
path_script = "R/"

# source functions ----
fili = dir(path_script,pattern="mw_[0123456789]")
lapply(X = fili, FUN = function(x) {source(paste0(path_script, x), echo=FALSE)})

# create savepoint repertory if not existing ----
dir.create(file.path("./", controls$savepoint), showWarnings = FALSE)

# common ----

# small custom functions ----

# aesthetics ----
