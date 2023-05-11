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
compute_cv = FALSE
# plot function
ppp_day = function(dat,mod,residuals=FALSE) {
  dd = mod$summary.fitted.values %>% 
    dplyr::bind_cols(dat) %>% 
    dplyr::mutate(res_mean=mean-reported_cases,
                  res_lwb=`0.025quant`-reported_cases,
                  res_upb=`0.975quant`-reported_cases)
  g1 =  dd %>% 
    ggplot(aes(x=date)) +
    geom_point(aes(y=reported_cases)) +
    geom_ribbon(aes(ymin=`0.025quant`,ymax=`0.975quant`),alpha=.5,fill=cust_cols[1]) +
    geom_line(aes(y=mean),colour=cust_cols[1]) +
    annotate(x=min(dat$date),y=min(dat$reported_cases),
             geom="text",
             label=paste0("WAIC=",round(mod$waic$waic)),hjust=0) 
  ll = max(abs(dd$res_mean))
  g2 =  dd %>% 
    ggplot(aes(x=date)) +
    geom_hline(yintercept=0,linetype=2) +
    geom_point(aes(y=res_mean),size=.7,shape=4,colour=cust_cols[1]) +
    coord_cartesian(ylim=c(-ll,ll)) +
    labs(y="Residuals")
  if(residuals) {
    g = cowplot::plot_grid(g1,g2,ncol=1,rel_heights = c(2,1))
  } else {
    g = g1
  }
  return(g)
}
# k-fold cross validation
kfoldcv = function(dat,mod,k=10,lk="nbinomial") {
  # extract formula
  ff = as.character(mod$.args$formula)
  ff = as.formula(paste(ff[2],ff[1],ff[3]))
  # adapt data
  dd = dat %>% 
    dplyr::mutate(k=sample(1:k,nrow(dat),replace=TRUE),
                  mean=NA,lbd=NA,upb=NA)
  # run k-fold cv
  for(i in 1:10) {
    rr = which(dd$k==i)
    ddx = dd %>% 
      dplyr::mutate(reported_cases = ifelse(k==i,NA,reported_cases))
    mm = INLA::inla(ff,
                    data = ddx,
                    family = lk,
                    control.mode = list(result = mod, restart=TRUE),
                    control.predictor = list(compute=TRUE,link=1))
    dd[rr,"mean"] = mm$summary.fitted.values[rr,"mean"]
    dd[rr,"lbd"] = mm$summary.fitted.values[rr,"0.025quant"]
    dd[rr,"upb"] = mm$summary.fitted.values[rr,"0.975quant"]
  }
  oo = dd %>% 
    dplyr::mutate(se=(reported_cases-mean)^2,
                  coverage=ifelse(reported_cases>lbd & reported_cases<upb,1,0)) %>% 
    summarise(rmse=sqrt(mean(se)),
              coverage=mean(coverage)) %>% 
    unlist() %>% unname()
  output = list(RMSE=round(oo[1],2),
                coverage=round(oo[2],2))
  return(output)
}
summary_exp = function(mod, pars) {
  mod$summary.fixed %>% 
    rownames_to_column() %>% 
    as_tibble() %>% 
    dplyr::filter(grepl(pars,rowname)) %>% 
    dplyr::transmute(rowname=rowname,
                     RR=exp(mean),
                     `0.025quant`=exp(`0.025quant`),
                     `0.975quant`=exp(`0.975quant`)) %>% 
    column_to_rownames()
}

#' We now attempt to link the measured SARS-CoV-2 viral load in wastewater to other indicators of SARS-CoV-2 infection in a population such as counts of laboratory-confirmed cases. Figure 1 shows the two time series and the scatter plot. A first issue is the question of the lag between these time series. Because of reporting delays, we would expect the reported cases to be delayed compared to the SARS-CoV-2 viral load in wastewater, although the duration of viral shedding may have an opposite effect. We use cross-correlation to assess the lag between reported cases and wastewater viral load (Figure 2). We find no clear indication of a time shift, and proceed by analyzing the comparative dynamics of both time series without a lag.
#' 
#+ fig.width=8, fig.height=5
# select one ARA
chosen_ara = "Aarau"
# data management
ww_one = ww1 %>% 
  filter(ara_name==chosen_ara) %>% 
  # remove missing
  filter(!is.na(vl),!is.na(reported_cases)) %>%
  # round reported_cases and hospitalizations
  mutate(reported_cases=round(reported_cases),
         reported_hospit=round(reported_hospit)) %>% 
  # log
  mutate(logvl=log2(vl)) %>% 
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
  geom_point(aes(x=date,y=reported_cases),colour=cust_cols[2]) +
  labs(x="Date",y="Reported cases")
g2b = ggplot(ww_one) +
  geom_point(aes(x=logvl,y=reported_cases),colour=cust_cols[2]) +
  # geom_smooth(aes(x=logvl,y=reported_cases),method="lm") +
  scale_y_continuous(trans="pseudo_log") +
  labs(y="Reported cases",x="Viral load (log)")
g2 = cowplot::plot_grid(g2a,g2b,rel_widths = c(2,1))
cowplot::plot_grid(g1,g2,labels=c("A","B"),ncol=1)

#' 
#' **Figure 1.** (A) SARS-CoV-2 viral load in wastewater in `r chosen_ara` over time and histogram. (B) Counts of laboratory-confirmed (or reported) cases of SARS-CoV-2 infection in the area covered by the  `r chosen_ara` wastewater plant and association with viral load. 
#'   
#+ fig.width=5, fig.height=3
tmp = cor.test(ww_one$reported_cases,ww_one$logvl,method = "pearson")
cc = tibble(lag=0,
            cor=tmp$estimate,
            cor_lwb=tmp$conf.int[1],
            cor_upb=tmp$conf.int[2])
for(i in 1:20) {
  dd = ww_one %>% 
    transmute(x=reported_cases,
              y=lead(logvl,i)) %>% 
    filter(!is.na(y))
  tmp = cor.test(dd$x,dd$y,method = "pearson")
  cc = cc %>% 
    bind_rows(tibble(lag=-i,
                     cor=tmp$estimate,
                     cor_lwb=tmp$conf.int[1],
                     cor_upb=tmp$conf.int[2]))
}
for(i in 1:20) {
  dd = ww_one %>% 
    transmute(x=reported_cases,
              y=lag(logvl,i)) %>% 
    filter(!is.na(y))
  tmp = cor.test(dd$x,dd$y,method = "pearson")
  cc = cc %>% 
    bind_rows(tibble(lag=i,
              cor=tmp$estimate,
              cor_lwb=tmp$conf.int[1],
              cor_upb=tmp$conf.int[2]))
}

cc %>% 
  ggplot() +
  geom_pointrange(aes(x=lag,y=cor,ymin=cor_lwb,ymax=cor_upb)) +
  labs(x="Lag of SARS-CoV-2 viral load compared to reported cases (days)",
       y="Pearson correlation (95% CI)") +
  annotate("segment",x=0.1,xend=3,y=0,yend=0,arrow=arrow(length=unit(0.1, "inches"))) +
  annotate("segment",x=-0.1,xend=-3,y=0,yend=0,arrow=arrow(length=unit(0.1, "inches"))) +
  annotate("text",x=2.5,y=.02,label="Reported cases behind",size=3,hjust=0) +
  annotate("text",x=-2.5,y=.02,label="Reported cases ahead",size=3,hjust=1)

#' 
#' **Figure 2.** Cross-correlation between SARS-CoV-2 reported cases and SARS-CoV-2 viral load in wastewater.
#'   

#' ## Model B1: multiplicative relation with laboratory-confirmed cases
#' 
#' We start with the simplest model, looking at a multiplicative association between viral load and reported cases. We use negative-binomial regression to allow for over-dispersion. This includes a log link, so that the number of reported cases on a given day $K_t$ is modelled as:
#' 
#' $$
#' \log(K_t) = \alpha + \beta log2V
#' $$
#' 
#' This 'log-log2' formulation is equivalent to a power law:
#' 
#' $$
#' K_t = \exp(\alpha) \times \exp(\beta)^{log2V}
#' $$
#' 
#'  We also include a week-end effect, as it has been shown that reported cases are lower during weekends.

mb1 = INLA::inla(reported_cases ~ 1 + 
                   weekend +
                   logvl,
                 data = ww_one,
                 family = "nbinomial",
                 control.compute = list(waic=TRUE))
summary(mb1)
summary_exp(mb1,pars="logvl|weekend")
ppp_day(ww_one,mb1)
if(compute_cv) kfoldcv(ww_one,mb1)

#' We find that indeed reported cases are about 35% lower during weekends. We also find that on average the number of reported cases on a given day is about 30% higher for every doubling in wastewater viral load on the same day. Still, the model fit is not very good: we observe 3 waves of decreasing amplitude over the period when looking at reported cases, while wastewater viral load suggests 4 waves of roughly similar amplitudes. This is likely explained by differences in ascertainment.
#' 
#' ## Model B2: multiplicative relation with laboratory-confirmed cases by period
#' 
#' We know that there has been changes in ascertainment over the period, with testing rates going down steadily during 2022, and even more after free-of-charge testing was lifted on January 1st, 2023. We thus decide to stratify the analysis using 4 period defined using the lowest counts of reported cases between waves and January 1st, 2023. Within these periods, we can reasonably consider ascertainment as stable, which allows us to analyse and compare dynamics. 
#' 
g1 = ggplot(ww_one) +
  geom_point(aes(x=date,y=logvl),colour=cust_cols[1]) +
  geom_vline(xintercept=cut_date, linetype=2) +
  labs(x="Date",y="Viral load (log)")
g2 = ggplot(ww_one) +
  geom_point(aes(x=date,y=reported_cases),colour=cust_cols[2]) +
  geom_vline(xintercept=cut_date, linetype=2) +
  labs(x="Date",y="Reported cases")
cowplot::plot_grid(g1,g2,ncol=1)
# ww_one %>% 
#   ggplot(aes(x=logvl,y=reported_cases)) +
#   geom_point() +
#   geom_smooth(method="lm") +
#   facet_grid(~period) +
#   scale_y_continuous(trans="pseudo_log")

mb2 = INLA::inla(reported_cases ~ 1 + 
                   weekend +
                   period +
                   period:logvl,
                 data = ww_one,
                 family = "nbinomial",
                 control.compute = list(waic=TRUE))
summary(mb2)
summary_exp(mb2,pars="logvl")
ppp_day(ww_one,mb2)
if(compute_cv) kfoldcv(ww_one,mb2)

#' 
#' We obtain a much better fit, with a clear association during the first period: a doubling of viral load corresponds to a multiplication of reported cases by 1.7, close to the value of 2 that would be expected if the relationship was perfect. For the other periods the RR decrease, arguably because of lower ascertainment. If we extrapolate, we can predict the counts of reported cases had ascertainment remained the same as during period 1, showing a much larger amplitude for waves 2, 3 and 4.

ww_one_p1 = ww_one %>% 
  dplyr::mutate(reported_cases=ifelse(period=="(-1,107]",reported_cases,NA))
mb2b = INLA::inla(reported_cases ~ 1 + 
                    weekend +
                    logvl,
                 data = ww_one_p1,
                 family = "nbinomial",
                 control.compute = list(waic=TRUE),
                 control.predictor=list(link=1))
ppp_day(ww_one,mb2b)

#' 
#' ## Model B3: non-linear relation
#' 
#' The log-linear
ww_one = ww_one %>% 
  dplyr::mutate(logvl_d=INLA::inla.group(logvl,n=15))
knots = cut_day
mb3 = INLA::inla(reported_cases ~ 1 + 
                   weekend +
                   ns(logvl, df = 10),
                 data = ww_one,
                 family = "nbinomial",
                 control.compute = list(waic=TRUE))
summary(mb3)
summary_exp(mb3,pars="logvl")
ppp_day(ww_one,mb3)
if(compute_cv) kfoldcv(ww_one,mb2)

