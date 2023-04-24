#' ---
#' title: "WBE for SARS-CoV-2 in Switzerland: model development (B) "
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
source("setup.R")
ww1 = readRDS(fs::path("../",controls$savepoint,"ww1.rds"))
# plot function
ppp_day = function(dat,mod) {
  mod$summary.fitted.values %>% 
    dplyr::bind_cols(dat) %>% 
    ggplot(aes(x=date)) +
    geom_point(aes(y=logvl)) +
    geom_ribbon(aes(ymin=`0.025quant`,ymax=`0.975quant`),alpha=.5,fill=cust_cols[1]) +
    geom_line(aes(y=mean),colour=cust_cols[1]) +
    annotate(x=min(dat$date),y=min(dat$logvl),geom="text",label=paste0("WAIC=",round(mod$waic$waic)),hjust=0)
}

#' We now attempt to link the measured SARS-CoV-2 viral load to other indicators of SARS-CoV-2 infection in a population such as counts of laboratory-confirmed cases and hospitalizations.
#' 
#+ fig.width=8, fig.height=7.5
# select one ARA
chosen_ara = "Aarau"
# data management
ww_one = ww1 %>% 
  filter(ara_name==chosen_ara) %>% 
  # remove missing
  filter(!is.na(vl)) %>%
  # log
  mutate(logvl=log(vl),
         logrep=log(reported_cases+1),
         loghos=log(reported_hospit+1)) %>% 
  # manage dates
  mutate(day=as.numeric(date)-min(as.numeric(date)),
         weekend=if_else(lubridate::wday(date,week_start=1)>=5,1,0)) %>% 
  rownames_to_column()
# plot
g1a = ggplot(ww_one) +
  geom_point(aes(x=date,y=logvl),colour=cust_cols[1]) +
  labs(x="Date",y="Viral load (log)")
g1b = ggplot(ww_one) +
  geom_histogram(aes(x=logvl),fill=cust_cols[1],colour="black",bins=30) +
  scale_y_continuous(expand=expansion(c(0,0.05))) +
  labs(x="Viral load (log)",y="Count")
g1 = cowplot::plot_grid(g1a,g1b,rel_widths = c(2,1))
g2a = ggplot(ww_one) +
  geom_point(aes(x=date,y=logrep),colour=cust_cols[2]) +
  labs(x="Date",y="Laboratory-confirmed cases (log)")
g2b = ggplot(ww_one) +
  geom_point(aes(x=logrep,y=logvl),colour=cust_cols[2]) +
  geom_smooth(aes(x=logrep,y=logvl),method="lm") +
  labs(x="Laboratory-confirmed cases (log)",y="Viral load (log)")
g2 = cowplot::plot_grid(g2a,g2b,rel_widths = c(2,1))
g3a = ggplot(ww_one) +
  geom_point(aes(x=date,y=loghos),colour=cust_cols[3]) +
  labs(x="Date",y="COVID-19 hospitalizations (log)")
g3b = ggplot(ww_one) +
  geom_point(aes(x=loghos,y=logvl),colour=cust_cols[3]) +
  geom_smooth(aes(x=loghos,y=logvl),method="lm") +
  labs(y="Viral load (log)",x="COVID-19 hospitalizations (log)")
g3 = cowplot::plot_grid(g3a,g3b,rel_widths = c(2,1))
cowplot::plot_grid(g1,g2,g3,labels=c("A","B","C"),ncol=1)

#' **Figure 1.** (A) SARS-CoV-2 viral load in wastewater in `r chosen_ara` over time and histogram. (B) Log counts of laboratory-confirmed cases of SARS-CoV-2 infection in the area covered by the  `r chosen_ara` wastewater plant and association with viral load. (C) Log counts of COVID-19 hospitalizations in the area covered by the  `r chosen_ara` wastewater plant and association with viral load.  
#'   
#' ## Model B1: log-linear relation with laboratory-confirmed cases
#' 
#' There seem to be a linear relation between log viral load and log laboratory-confirmed cases.

mb1 = INLA::inla(logvl ~ 1 + 
                   logrep,
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE))
summary(mb1)
ppp_day(ww_one,mb1)

# 
# # model A7: reported cases by period (defined using low points in reported cases)
# cut_date = ymd("2022-05-15","2022-09-15","2023-01-01")
# cut_day = as.numeric(cut_date) - min(as.numeric(ww_one$date))
# ww_one = ww_one %>% 
#   dplyr::mutate(period=cut(day,breaks=c(-1,cut_day,max(ww_one$day)+1)))
# cowplot::plot_grid(ggplot(ww_one) + geom_point(aes(x=date,y=logrep),col="red") + geom_vline(xintercept=cut_date,linetype=2),
#                    ggplot(ww_one) + geom_point(aes(x=date,y=logvl)) + geom_vline(xintercept=cut_date,linetype=2),
#                    ncol=1)
# ggplot(ww_one) + geom_point(aes(x=logrep,y=logvl,col=period)) + facet_wrap(~period)
# ma7 = INLA::inla(logvl ~ 1 + 
#                    logrep:period, 
#                  data = ww_one, 
#                  family = "gamma",
#                  control.compute = list(waic=TRUE))
# summary(ma7)
# ppp_day(ww_one,ma7) + geom_vline(xintercept=cut_date,linetype=2)
# 
# ww_one %>% 
#   bind_cols(ma7$summary.fitted.values) %>% 
#   ggplot(aes(x=logrep)) + 
#   geom_point(aes(y=logvl,col=period)) + 
#   geom_ribbon(aes(ymin=`0.025quant`,ymax=`0.975quant`),alpha=.5,fill=cust_cols[1]) +
#   geom_line(aes(y=mean),colour=cust_cols[1]) +
#   facet_wrap(~period)
# 
# 
# 
# ww_one %>% 
#   mutate(mean=predict(ma7b)) %>% 
#   ggplot(aes(x=logrep)) + 
#   geom_point(aes(y=logvl,col=period)) + 
#   # geom_ribbon(aes(ymin=`0.025quant`,ymax=`0.975quant`),alpha=.5,fill=cust_cols[1]) +
#   geom_line(aes(y=mean),colour=cust_cols[1]) +
#   facet_wrap(~period)
# 
# 
# # model A8: reported cases by period + remaining time trend
# ma8 = INLA::inla(logvl ~ 1 + 
#                    logrep:period +
#                    f(day,period,
#                      model="rw2",
#                      hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.001)))), 
#                  data = ww_one, 
#                  family = "gaussian",
#                  control.compute = list(waic=TRUE))
# summary(ma8)
# ppp_day(ww_one,ma8)
# 
