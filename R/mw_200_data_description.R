#' ---
#' title: "WBE for SARS-CoV-2 in Switzerland: data description "
#' author: "Julien Riou"
#' date: "`r Sys.Date()`"
#' params:
#'    controls: controls
#' output:
#'    html_document:
#'      code_folding : hide
#' toc: true
#' toc_float: true
#' toc_depth: 4
#' number_sections: true
#' highlight: pygments
#' theme: cosmo
#' link-citations: true
#' ---


#+ results="hide", warnings="false", echo="false"
# run setup
source("setup.R")
# load data
ww1 = readRDS(fs::path("../",controls$savepoint,"ww1.rds"))

#' We use measurements of SARS-CoV-2 concentration in wastewater from multiple ARAs in Switzerland in 2022 and 2023. Viral concentration (*C*, unit: gene copies [gc] per liter) is transformed into viral load (*V*; unit: gc per day per 100,000) using the flow of wastewater on the same day (*F*) and the size of the population covered (*P*):
#' $$
#'  V = \frac{C \times F}{P/100,000}
#' $$
#' Table 1 provides a summary of the available data. 
#' Reporting frequencies and periods depended on the ARA, with daily measurements for the full period only available in a few cases (Figure 1). ARAs also sent their samples to different laboratories. In some cases, there were also changes in the method used.
#' SARS-CoV-2 could be detected in the wastewater in most cases, with a few occurrences of no detection (Figure 2).
#' Viral load varied over time, with large heterogeneity across ARAs, although some patterns emerge on visual inspection (Figure 3). 
#' 
#' **Table 1.** Summary of available data.
mw_100_desc_table(ww1) %>% 
  dplyr::mutate(across(everything(),as.character)) %>% 
  tidyr::gather() %>% 
  flextable::flextable(cwidth=c(3,4)) 

#+ fig.width=8, fig.height=10
mw_101_fig_missing(ww1)
#' **Figure 1.** Available measurements over time by ARA (grouped by canton).

#+ fig.width=8, fig.height=10
mw_103_fig_detect(ww1,detect_limit=0)
#' **Figure 2.** SARS-CoV-2 detection in wastewater over time by ARA (limit=0).

#+ fig.width=8, fig.height=10
mw_104_fig_vl(ww1,lower_limit=5000)
#' **Figure 3.** Weekly mean SARS-CoV-2 viral load in wastewater by ARA (removing one outlier below 5,000).

