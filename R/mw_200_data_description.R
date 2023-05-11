#' ---
#' title: "WBE for SARS-CoV-2 in Switzerland: data description "
#' author: "Julien Riou"
#' date: "`r Sys.Date()`"
#' params:
#'    controls: controls
#' output:
#'    html_document:
#'      code_folding : hide
<<<<<<< HEAD
#' toc: true
#' toc_float: true
#' toc_depth: 4
#' number_sections: true
#' highlight: pygments
#' theme: cosmo
#' link-citations: true
=======
#'      toc: true
#'      toc_float: true
#'      toc_depth: 4
#'      number_sections: true
#'      highlight: pygments
#'      theme: cosmo
#'      link-citations: true
>>>>>>> gh-pages
#' ---


#+ results="hide", warnings="false", echo="false"
# run setup
source("setup.R")
# load data
ww1 = readRDS(fs::path("../",controls$savepoint,"ww1.rds"))
<<<<<<< HEAD
=======
shapes = readRDS(fs::path("../",controls$savepoint,"shapes.rds"))
>>>>>>> gh-pages

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

<<<<<<< HEAD
#+ fig.width=8, fig.height=10
mw_103_fig_detect(ww1,detect_limit=0)
#' **Figure 2.** SARS-CoV-2 detection in wastewater over time by ARA (limit=0).

#+ fig.width=8, fig.height=10
mw_104_fig_vl(ww1,lower_limit=5000)
#' **Figure 3.** Weekly mean SARS-CoV-2 viral load in wastewater by ARA (removing one outlier below 5,000).
=======
#+ fig.width=8, fig.height=6
mw_110_map_missing(ww1,shapes)
#' **Figure 2.** Total measurements by ARA.

#+ fig.width=8, fig.height=10
mw_103_fig_detect(ww1,detect_limit=0) # TODO: replace by something based on LOQ/LOD
#' **Figure 3.** SARS-CoV-2 detection in wastewater over time by ARA.

#+ fig.width=8, fig.height=10
mw_104_fig_vl(ww1,lower_limit=5000)
#' **Figure 4.** Weekly mean SARS-CoV-2 viral load in wastewater by ARA (removing one outlier below 5,000). Dashed lines show the delimitation in four periods.

#+ fig.width=8, fig.height=6
mw_111_map_vl(ww1,shapes,lower_limit=5000)
#' **Figure 5.** Mean SARS-CoV-2 viral load in wastewater by ARA by period.



if(FALSE) {
  ggplot(ww1) +
    geom_point(aes(x=pop,y=pop_total)) +
    geom_abline(intercept=0,slope=1)
  
  dd = ww1 %>% 
    group_by(ara_id,ara_name) %>% 
    summarise(pop=max(pop),
              pop_total=max(pop_total)) %>% 
    mutate(rel=(pop_total-pop)/pop,
           abs=pop-pop_total) %>% 
    arrange(-abs(rel))
  ggplot(dd) +
    geom_point(aes(x=ara_name,y=pop),colour="firebrick")+
    geom_point(aes(x=ara_name,y=pop_total),colour="dodgerblue",shape=1) +
    scale_x_discrete(limits=rev) +
    coord_flip() +
  theme(axis.text = element_text(size=5))
}
>>>>>>> gh-pages

