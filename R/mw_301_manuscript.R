#' ---
#' title: "Tables and figures for manuscript"
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
labo_names = tibble(lab_method=c("EAWAG_0","ALTGR_0","Eurofins_0","Eurofins_1","KLZH_0","LdU_0","SCAV_NE_0","Microsynth_0","Microsynth_1","SUPSI_0"),
                 labo_label=c("A (ref.)","B","C1","C2","D","E","F","G1","G2","H"))
labo_cov = tibble(Variable=c("lab_methodALTGR_0","lab_methodEurofins_0","lab_methodEurofins_1","lab_methodKLZH_0",      
                            "lab_methodLdU_0","lab_methodSCAV_NE_0","lab_methodMicrosynth_0" ,"lab_methodMicrosynth_1",  
                            "lab_methodSUPSI_0","weekend","hol","prop_under_20b",        
                            "prop_over_65b","ssep3_medb","employment_factorb"),
                  labo_label=c("B","C1","C2","D","E","F","G1","G2","H",
                       "Week-end","Public holiday","Proportion aged <20","Proportion aged >65","Median Swiss-SEP","Employment factor"))
ww_all = ww_all %>% 
  left_join(wwtpnumber,by="ara_kt") %>% 
  left_join(labo_names,by="lab_method") %>%   
  dplyr::mutate(ara_nuts=paste0(NUTS2_name," / ",ara_name)) %>% 
  left_join(wwtpnutsnumber,by="ara_nuts")


# Description -------------------------------------------------------------

mw_100_desc_table(ww_all) %>%
  dplyr::mutate(across(everything(),as.character)) %>%
  tidyr::gather() %>%
  xtable::xtable() %>% 
  print(include.rownames = FALSE)

ww_all %>% 
  group_by(ara_name) %>% 
  summarise(n=n(),
            min=min(date),
            max=max(date),
            dur=as.numeric(max-min),
            freq=n/dur*7) %>% view()

ww_all %>% 
  filter(below_lod==0,below_loq==0) %>% 
  summarise(n=n(),
            min=formatC(min(vl),format="e",digits=1),
            max=formatC(max(vl),format="e",digits=1),
            sd=var(vl))
7.4e+14 /  1.7e+10

x = ww_all %>% 
  group_by(wwtp_index) %>% 
  summarise(m=mean(vl)) %>% 
  summarise(min=min(m),
            max=max(m))
x$max/x$min


ww_all$kt %>% unique() %>% length()

ww_all %>% 
  group_by(period) %>% 
  summarise(nara=length(unique(ara_id)))

ww_all %>% 
  group_by(kt) %>% 
  summarise(min=min(wwtp_index),
            max=max(wwtp_index)) %>% 
  print(n=26)


# Modèle ------------------------------------------------------------------

source("inla.rsquared.R")
inla.rsquared(ww_all,ma5.4.5)
summary(ma5.4.5)

# Figure 1 ------------------------------------------------------------------------------

g1 = mw_110_map_missing(ww_all,shapes)
g2 = mw_106_fig_vl_time(ww_all)
g3 = mw_104_fig_vl_nuts(ww_all)

leftp = cowplot::plot_grid(g1,g2,labels=c("A","B"),ncol=1,rel_heights = c(1.5,1))
f1 = cowplot::plot_grid(leftp,g3,labels=c("","C"))
ggsave(plot=f1, file=fs::path("../",controls$savepoint,"fig1.pdf"),width=9,height=6.3)

# Figure 2 ------------------------------------------------------------------------------

summary_exp_vl(ma5.4.5,pars="lab|method|hol|weekend|pop_total|prop_under_20|prop_over_65|prop_non_ch_eu|ssep3_|emp")
g1 = plot_exp_vl(ma5.4.5,pars="lab", ref="A (ref.)",labs=labo_cov,col="firebrick")+ coord_cartesian(ylim=c(0,2))
g2 = plot_exp_vl(ma5.4.5,pars="hol|weekend",labs=labo_cov,col="dodgerblue") + coord_cartesian(ylim=c(.7,1.3))
g3 = plot_exp_vl(ma5.4.5,pars="pop_total|prop_under_20|prop_over_65|prop_non_ch_eu|ssep3_|emp",labs=labo_cov,col="seagreen")   + coord_cartesian(ylim=c(.7,1.3))

f2 = cowplot::plot_grid(g1,g3,g2,labels=c("A","B","C"),align="hv",ncol=3,rel_widths = c(2,2,1.2))
ggsave(plot=f2, file=fs::path("../",controls$savepoint,"fig2.pdf"),width=9,height=3.5)


ppp_vl_ara(ww_all,ma5.4.5)
#+ ma5.4.1b, fig.width=8, fig.height=4,  R.options = list(width = 1000)
avg_time_trend_reg(ww_all,ma5.4.5)
#+ ma5.4.1c, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_130_map_relative_vl(ma5.4.5,corr_all_ara,shapes)
#+ ma5.4.1d, fig.width=8, fig.height=6,  R.options = list(width = 1000)
mw_131_map_deviation_from_average(ma5.4.5,corr_all_ara,ww_all,shapes,12)


# Figure 3 ------------------------------------------------------------------------------


# region_labs = tibble(nam=c("Northwest","Central","Lake Geneva","Zurich","Eastern","Mittelland","Ticino"),
#                      x=c("2022-05-30","2022-07-10","2023-02-10","2023-04-01","2023-09-20","2023-09-10","2023-09-30"),
#                      y=c(10,.14,.16,5,8,5,.09),
#                      xend=c("2022-03-14","2022-05-24","2023-01-26","2023-03-22","2023-11-05","2023-11-03","2023-08-13"),
#                      yend=c(12,.12,.24,3.6,6.5,4.7,.11))
region_labs = tibble(nam=c("Northwest","Northwest","Central","Zurich","Eastern","Mittelland"),
                     x=c("2022-05-25","2022-05-25","2022-07-10","2023-04-01","2023-09-20","2023-09-10"),
                     y=c(9,9,.14,5.5,8,5),
                     xend=c("2022-03-14","2022-07-17","2022-05-24","2023-03-22","2023-11-12","2023-11-08"),
                     yend=c(12,9,.13,3.5,5.5,2.9))

g1 = avg_time_trend_reg(ww_all,ma5.4.5) +
  geom_segment(data=region_labs,aes(x=as.Date(x),y=y,xend=as.Date(xend),yend=yend),linewidth=.2) +
  geom_label(data=region_labs,aes(x=as.Date(x),y=y,label=nam,colour=nam),size=3,show.legend=FALSE)

g2 = mw_140_regional_correlation(ww_all,ma5.4.5,corr_all_ara,type="hospitalizations") 

g3 = mw_130_map_relative_vl(ma5.4.5,corr_all_ara,shapes)

botpn = cowplot::plot_grid(g2,g3,labels=c("B","C"),align="v",ncol=2,rel_widths = c(1,1.6))
f3 = cowplot::plot_grid(g1,botpn,ncol=1,labels=c("A",""),rel_heights = c(1,1.6))
ggsave(plot=f3, file=fs::path("../",controls$savepoint,"fig3.pdf"),width=9,height=7.3)


mw_130_map_relative_vl(ma5.4.5,corr_all_ara,shapes,pprint=TRUE)

supf = mw_142_regional_correlation_lag(dat=ww_all,
                                mod=ma5.4.5,
                                corr=corr_all_ara,
                                type="hospitalizations") 
mw_143_regional_hospits(ww_all)
