#' ---
#' title: "WBE for SARS-CoV-2 in Switzerland: model development (A) "
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

#' We start model development focusing on one ARA, *Aarau*. We consider measurements of SARS-CoV-2 viral load on the logarithmic scale. The objective of this first part is to create a model able to describe the variation in SARS-CoV-2 viral load in wastewater over time in one location, depending on covariates, that can be later extended to multiple location. We use a Bayesian approach with *integrated nested Laplace approximation* as implemented in the `R-INLA` package, as it integrates many tools for spatial modelling. For each fitted model we present the model output, posterior predictive plots and a measure of the goodness-of-fit using the WAIC (Watanabe–Akaike information criterion).

#+ fig.width=8, fig.height=3.5
# select one ARA
chosen_ara = "Aarau"
# data management
ww_one = ww %>% 
  filter(ara_name==chosen_ara) %>% 
  # remove missing
  filter(!is.na(vl), vl>lower_limit) %>%
  # log
  mutate(logvl=log(vl),
         logrep=log(reported_cases+1)) %>% 
  # manage dates
  mutate(day=as.numeric(date)-min(as.numeric(date))) %>% 
  rownames_to_column()
# plot
g1 = ggplot(ww_one) +
  geom_point(aes(x=date,y=logvl),colour=cust_cols[1]) +
  labs(x="Date",y="Viral load (log)")
g2 = ggplot(ww_one) +
  geom_histogram(aes(x=logvl),fill=cust_cols[1],colour="black",bins=30) +
  scale_y_continuous(expand=expansion(c(0,0.05))) +
  labs(x="Viral load (log)",y="Count")
cowplot::plot_grid(g1,g2,labels=c("A","B"),rel_widths = c(2,1))
#' **Figure 1.** SARS-CoV-2 viral load in wastewater in `r chosen_ara` (panel A: measurements over time; panel B: histogram).

#'   
#' ## Model A1: gamma regression
#' 
#' We select gamma regression with a log link as the dependent variable is continuous, strictly positive and on the log scale. The log link also implies that all effects will be multiplicative.
#' 
#' The first, basic model only includes a linear trend over time.
ma1 = INLA::inla(logvl ~ 1 + day,
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE))
summary(ma1)
#+ fig.width=4, fig.height=3
ppp_day(ww_one,ma1)

# 
# # model A2: splines
# dknots = seq(from=0,to=max(ww_one$day),by=30)
# ma2 = INLA::inla(logvl ~ 1 + 
#                    bs(day,knots=dknots), 
#                  data = ww_one, 
#                  family = "gamma",
#                  control.compute = list(waic=TRUE))
# summary(ma2)
# ppp_day(ww_one,ma2)
# ma2$waic$waic
# 
# # model A3: rw1
# ma3 = INLA::inla(logvl ~ 1 + 
#                    f(day,
#                      model="rw1",
#                      hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))), 
#                  data = ww_one, 
#                  family = "gamma",
#                  control.compute = list(waic=TRUE))
# summary(ma3)
# ppp_day(ww_one,ma3)
# 
# # model A4: rw2
# ma4 = INLA::inla(logvl ~ 1 + 
#                    f(day,
#                      model="rw2",
#                      hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))), 
#                  data = ww_one, 
#                  family = "gamma",
#                  control.compute = list(waic=TRUE))
# summary(ma4)
# ppp_day(ww_one,ma4)
# 
# # model A5: rw2 + method change
# ma5 = INLA::inla(logvl ~ 1 + 
#                    f(day,
#                      model="rw2",
#                      hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#                    method, 
#                  data = ww_one, 
#                  family = "gamma",
#                  control.compute = list(waic=TRUE))
# summary(ma5)
# ppp_day(ww_one,ma5)
# 
# # model A6: reported cases
# cowplot::plot_grid(ggplot(ww_one) + geom_point(aes(x=date,y=logrep)),
#                    ggplot(ww_one) + geom_point(aes(x=date,y=logvl),col="red"),
#                    ncol=1)
# ma6 = INLA::inla(logvl ~ 1 + logrep, 
#                  data = ww_one, 
#                  family = "gamma",
#                  control.compute = list(waic=TRUE))
# summary(ma6)
# ppp_day(ww_one,ma6)
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
