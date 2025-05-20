#' ---
#' title: "S1 text"
#' author: "Julien Riou"
#' date: "`r Sys.Date()`"
#' ---



# Select ------------------------------------------------------------------
controls = readRDS(file.path("savepoints/savepoint_2025-01-24/controls.rds"))
source("R/setup.R")
setwd("R")
shapes = readRDS(fs::path("../",controls$savepoint,"shapes.rds"))
ww_all = readRDS(file=paste0("../",controls$savepoint,"ww_all.rds"))
savepath = "C:/Users/ju6558/Documents/en_cours/2024_wbe/667142f8f7eb62e0d8010282"
ma5.4.5 = readRDS(fs::path("../",controls$savepoint,"ma5.4.5.rds"))
corr_all_ara = readRDS(file=paste0("../",controls$savepoint,"corr_all_ara.rds"))


# Relabelling -------------------------------------------------------------

# WWTP number from 1 to 118 sorted by canton
wwtpnumber = ww_all %>% 
  dplyr::filter(below_lod==0,below_loq==0) %>% 
  dplyr::group_by(kt,ara_kt,week) %>% 
  dplyr::summarise(mvl=median(vl,na.rm=TRUE),.groups="drop") %>% 
  dplyr::mutate(wwtp_index=as.numeric(as.factor(ara_kt))) %>% 
  dplyr::select(ara_kt,wwtp_index) %>% 
  dplyr::distinct() 
wwtpnutsnumber = ww_all %>% 
  dplyr::filter(below_lod==0,below_loq==0) %>% 
  dplyr::mutate(ara_nuts=paste0(NUTS2_name," / ",ara_name)) %>% 
  dplyr::group_by(NUTS2_name,kt,ara_nuts,week) %>% 
  dplyr::summarise(mvl=median(vl,na.rm=TRUE),.groups="drop") %>% 
  dplyr::mutate(wwtp_nuts_index=as.numeric(as.factor(ara_nuts))) %>% 
  dplyr::select(ara_nuts,wwtp_nuts_index) %>% 
  dplyr::distinct() 


# lab numbers as letters with EAWAG first
labo_names = tibble(lab_method=c("EAWAG_0","ALTGR_0","Eurofins_0","Eurofins_1","KLZH_0","LdU_0","Microsynth_0","Microsynth_1","SCAV_NE_0","SUPSI_0"),
                    labo_label=c("A (ref.)","B","C1","C2","D","E","F","G1","G2","H"))
labo_cov = tibble(Variable=c("lab_methodALTGR_0","lab_methodEurofins_0","lab_methodEurofins_1","lab_methodKLZH_0",      
                             "lab_methodLdU_0","lab_methodMicrosynth_0" ,"lab_methodMicrosynth_1", "lab_methodSCAV_NE_0",   
                             "lab_methodSUPSI_0","weekend","hol","prop_under_20b",        
                             "prop_over_65b","ssep3_medb","employment_factorb"),
                  labo_label=c("B","C1","C2","D","E","F","G1","G2","H",
                               "Week-end","Public holiday","Proportion aged <20","Proportion aged >65","Median Swiss-SEP","Employment factor"))
ww_all = ww_all %>% 
  left_join(wwtpnumber,by="ara_kt") %>% 
  left_join(labo_names,by="lab_method") %>%   
  dplyr::mutate(ara_nuts=paste0(NUTS2_name," / ",ara_name)) %>% 
  left_join(wwtpnutsnumber,by="ara_nuts")

#' # Supplementary methods

# methods model

# model fit

#' 
#' # Correlation with hospitalisations
#' 

g1 = mw_141_crude_correlation(ww_all,corr_all_ara,pprint=FALSE)
g2 = mw_140_regional_correlation(ww_all,ma5.4.5,corr_all_ara,type="hospitalizations",pprint=FALSE) 

cowplot::plot_grid(g1,g2,align = "hv",labels = LETTERS)
ggsave(file="../supplementary/correlation_hospitalisation.png",width = 8, height=4)

corr1 = mw_141_crude_correlation(ww_all,corr_all_ara,pprint=TRUE)
corr2 = mw_140_regional_correlation(ww_all,ma5.4.5,corr_all_ara,type="hospitalizations",pprint=TRUE) 

corr = corr1 %>% 
  rename(cor1=cor) %>% 
  left_join(corr2,by = join_by(NUTS2_name, period)) %>% 
  rename(cor2=cor) %>% 
  ungroup()

corr %>% 
  mutate(inc=cor2>cor1) %>% 
  summarise(increase=mean(inc))

corr %>% 
  summarise(cor1=median(cor1),cor2=median(cor2))

corr %>% 
  group_by(period) %>% 
  summarise(cor1=median(cor1),cor2=median(cor2))


ww_all %>%
  group_by(period,lab_method) %>% 
  count() %>% view()


mw_142_regional_correlation_lag(dat=ww_all,
                                       mod=ma5.4.5,
                                       corr=corr_all_ara,
                                       type="hospitalizations") 
ggsave(file="../supplementary/correlation_lag.png",width = 8, height=4)

# coeff of determination


#' 
#' 
#' #' We use measurements of SARS-CoV-2 concentration in wastewater from multiple ARAs in Switzerland in 2022 and 2023. 
#' #' Viral concentration (*C*, unit: gene copies [gc] per liter) is transformed into viral load (*V*; unit: gc per day per 100,000) 
#' #' using the flow of wastewater on the same day (*F*) and the size of the population covered (*P*):
#' #' $$
#' #'  V = \frac{C \times F}{P/100,000}
#' #' $$
#' #' Table 1 provides a summary of the available data. 
#' #' In total 118 ARAs reported data over the full period.
#' #' Reporting frequencies and periods varied depending on the ARA, with daily measurements for the full period only available 
#' #' in a few cases (Figure 1). ARAs also sent their samples to different laboratories. In some cases, there were also changes 
#' #' in the laboratory method used.
#' #' SARS-CoV-2 could be detected in the wastewater in most cases, with a few occurrences of no detection (below LOD) and some measurements
#' #' below the limit of quantification (LOQ), leading to higher measurement error (Figure 2).
#' #' Viral load varied over time, with large heterogeneity across ARAs, although some patterns emerge on visual inspection (Figure 3). 
#' #' We also consider covariates at the ARA level. While population is already accounted for in the computation of 
#' #' the viral load, covariates included total population covered by the ARA. We also considered the proportion aged below 20 
#' #' or above 65, the proportion of non-Swiss non-EU citizen, the average Swiss index of socio-economic position (weighted by hectare),
#' #' and the *employment factor* (working places divided by population).
#' #' 
#' 
#' #+ desc_tab
#' mw_100_desc_table(ww1) %>% 
#'   dplyr::mutate(across(everything(),as.character)) %>% 
#'   tidyr::gather() %>% 
#'   dplyr::rename(Variable=key,Value=value) %>% 
#'   flextable::flextable(cwidth=c(4,10)) 
#' #' **Table 1.** Summary of available data.
#' 
#' #+ desc_tab_nuts2
#' mw_100_desc_table(ww1,NUTS2_name) %>% 
#'   dplyr::select(1:7)  %>% 
#'   flextable::flextable(cwidth=rep(4,7)) 
#' #' **Table 2.** Summary of available data by NUTS-2 region.
#' 
#' #+ fig_missing, fig.width=8, fig.height=10
#' mw_101_fig_missing(ww1)
#' #' **Figure 1.** Available measurements over time by ARA (grouped by canton).
#' 
#' #+ map_missing, fig.width=8, fig.height=6
#' mw_110_map_missing(ww1,shapes)
#' #' **Figure 2.** Total measurements by ARA.
#' 
#' #+ fig_detect, fig.width=8, fig.height=10
#' mw_103_fig_detect(ww1)
#' #' **Figure 3.** SARS-CoV-2 detection in wastewater over time by ARA.
#' 
#' #+ fig_vl2, fig.width=8, fig.height=6
#' mw_106_fig_vl_time(ww1)
#' #' **Figure 4.** Daily SARS-CoV-2 viral load in wastewater by ARA (removing values below the LOD or LOQ).
#' 
#' #+ fig_vl, fig.width=8, fig.height=10
#' mw_104_fig_vl(ww1)
#' #' **Figure 5.** Weekly median SARS-CoV-2 viral load in wastewater by ARA (removing values below the LOD or LOQ). Dashed lines show the delimitation in four periods.
#' 
#' #+ map_vl, fig.width=8, fig.height=14
#' mw_111_map_vl(ww1,shapes)
#' #' **Figure 6.** Median SARS-CoV-2 viral load in wastewater by ARA by period.
#' 
#' 
#' 
#' 
