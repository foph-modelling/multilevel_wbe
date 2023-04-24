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
# k-fold cross validation
kfoldcv = function(dat,mod,k=10) {
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
      dplyr::mutate(logvl = ifelse(k==i,NA,logvl))
    mm = INLA::inla(ff,
               data = ddx,
               family = "gamma",
               control.mode = list(result = mod, restart=TRUE),
               control.predictor = list(compute=TRUE,link=1))
    dd[rr,"mean"] = mm$summary.fitted.values[rr,"mean"]
    dd[rr,"lbd"] = mm$summary.fitted.values[rr,"0.025quant"]
    dd[rr,"upb"] = mm$summary.fitted.values[rr,"0.975quant"]
  }
  oo = dd %>% 
    dplyr::mutate(se=(logvl-mean)^2,
                  coverage=ifelse(logvl>lbd & logvl<upb,1,0)) %>% 
    summarise(rmse=sqrt(mean(se)),
              coverage=mean(coverage)) %>% 
    unlist() %>% unname()
  output = list(RMSE=round(oo[1],2),
                coverage=round(oo[2],2))
  return(output)
}

#' We start model development focusing on one ARA. We consider measurements of SARS-CoV-2 viral load on the logarithmic scale. The objective of this first part is to create a model able to describe the variation in SARS-CoV-2 viral load in wastewater over time in one location, depending on covariates, that can be later extended to multiple location. We use a Bayesian approach with *integrated nested Laplace approximation* as implemented in the `R-INLA` package, as it integrates many tools for spatial modelling. For each fitted model we present the model output, posterior predictive plots and a measure of the goodness-of-fit using the WAIC (Watanabe–Akaike information criterion). We also conduct 10-fold cross-validation and consider the root-mean squared error and the coverage in the cross-validation as the main indicator for model selection.

#+ fig.width=8, fig.height=3.5
# select one ARA
chosen_ara = "Aarau"
# data management
ww_one = ww1 %>% 
  filter(ara_name==chosen_ara) %>% 
  # remove missing
  filter(!is.na(vl)) %>%
  # log
  mutate(logvl=log(vl)) %>% 
  # manage dates
  mutate(day=as.numeric(date)-min(as.numeric(date)),
         weekend=if_else(lubridate::wday(date,week_start=1)>=5,1,0)) %>% 
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
#' We select gamma regression with a log link as the dependent variable is continuous and strictly positive. The log link also implies that all effects will be multiplicative. This is espacially useful as ultimately changes in the log viral load over time are expected to follow the epidemic dynamic of SARS-CoV-2, so with exponential growth and exponential decay. 
#' 
#' The first, basic model only includes a linear trend over time. 
#' 
#' $$
#' \log(logV) = \alpha + \beta t
#' $$
#' 

ma1 = INLA::inla(logvl ~ 1 + day,
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE,config=TRUE))
summary(ma1)
kfoldcv(ww_one,ma1)
#+ fig.width=4, fig.height=3, fig.align="center"
ppp_day(ww_one,ma1)

#' In the output we can see the call, the used time in seconds, the posterior estimates of the parameters and hyperparameters, and the WAIC. We also look at the posterior predictive plot comparing the data with model predictions. Obviously it will not produce a good fit of the data, it's just meant as a starting point. 
#' 
#' 
#'   
#' ## Model A2: random walk
#' 
#' We deal with the variation in log viral load over time with a random walk, so that the difference between two successive observations follows a normal distribution.

ma2 = INLA::inla(logvl ~ 1 +
                   f(day, model="rw1", scale.model=TRUE, constr=TRUE),
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE))
summary(ma2)
kfoldcv(ww_one,ma2)
ppp_day(ww_one,ma2)

#' We see that the model fits much better, although visual inspection indicates over-fitting.
#' 
#' ## Model A3: random walk of order 2
#' 
#' We attempt to reduce over-fitting by using a random-walk of order 2, so that the difference depends on the last 2 observations.

ma3 = INLA::inla(logvl ~ 1 +
                   f(day,model="rw2", scale.model=TRUE, constr=TRUE),
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE))
summary(ma3)
kfoldcv(ww_one,ma3)
ppp_day(ww_one,ma3)

#' We still observe some over-fitting, even though there is a clear improvement.
#' 
#' ## Model A4: random walk of order 2 + penalized complexity
#' 
#' We reduce over-fitting by adding priors that penalize complexity (Simpson et al, 2017).
#'  
ma4 = INLA::inla(logvl ~ 1 +
                   f(day,model="rw2", scale.model=TRUE, constr=TRUE,
                     hyper=list(prec = list(prior = "pc.prec", param = c(0.1, 0.01)))),
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE))
summary(ma4)
kfoldcv(ww_one,ma4)
ppp_day(ww_one,ma4)

#' The model now appropriately captures the dynamics of SARS-CoV-2 viral load, although still underestimates the variability especially during certain periods.
#'  
#' ## Model A5: additional covariates
#'  
#' There are other covariates available at the level of an ARA that may explain part of the residual variability. One is a change in the methods used to quantify SARS-CoV-2. A second is the day of the week, as behaviour may change for instance during weekends.

ggplot(ww_one) +
  geom_point(aes(x=date,y=logvl)) +
  geom_line(aes(x=date,y=31,colour=factor(method))) +
  labs(colour="Method")

ma5 = INLA::inla(logvl ~ 1 +
                   f(day,model="rw2", scale.model=TRUE, constr=TRUE,
                     hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                   method +
                   weekend,
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE))
summary(ma5)
kfoldcv(ww_one,ma5)
ppp_day(ww_one,ma5)

#' In this case at least, we don't see much of an impact of these two covariates.
#' 
#' ## Model A6: AR
#'  
#' A similar idea is to use an auto-regressive model of order 1, where the value of observation $i$ depends on the value of observation $i-1$ scaled by a parameter $\rho$.
#' 
ma6 = INLA::inla(logvl ~ 1 +
                   f(day,model="ar1", constr=TRUE),
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE))
summary(ma6)
kfoldcv(ww_one,ma6)
ppp_day(ww_one,ma6)

#' Results appear similar to model A2.
#' 
#' ## Model A7: splines
#'  
#' An alternative could be using splines. The downside is that we have to choose knots. We can start with a knot every 30 days.
#' 
cknots = seq(from=0,to=max(ww_one$day),by=30)
ma7 = INLA::inla(logvl ~ 1 +
                   bs(day,knots = cknots),
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE))
summary(ma7)
kfoldcv(ww_one,ma7)
ppp_day(ww_one,ma7)

#' Results appear similar to model A4.
#' 

