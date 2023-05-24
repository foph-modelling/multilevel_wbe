#' ---
#' title: "WBE for SARS-CoV-2 in Switzerland: model development (A1) "
#' author: "Julien Riou"
#' date: "`r Sys.Date()`"
#' params:
#'    controls: controls
#' output:
#'    html_document:
#'      code_folding : hide
#'      toc: true
#'      toc_float: true
#'      toc_depth: 4
#'      number_sections: false
#'      highlight: pygments
#'      theme: cosmo
#'      link-citations: true
#' ---


#+ results="hide", warnings="false", echo="false"
source("setup.R")
ww1 = readRDS(fs::path("../",controls$savepoint,"ww1.rds"))
shapes = readRDS(fs::path("../",controls$savepoint,"ww1.rds"))

#' # One ARA
#' 
#' We start model development focusing on one ARA. The objective of this first part is to create a model able to describe the variation in SARS-CoV-2 viral load in wastewater over time in one location depending on covariates. The model will later be extended to multiple location. We use a Bayesian approach based on *integrated nested Laplace approximation* as implemented in the `R-INLA` package, as it integrates many tools for spatial modelling. For each fitted model we present the model output, posterior predictive plots and a measure of the goodness-of-fit using the WAIC (Watanabe–Akaike information criterion). We also conduct 10-fold cross-validation and consider the root-mean squared error and the coverage in the cross-validation as the main indicator for model selection.

#+ fig.width=8, fig.height=3.5
# select one ARA
chosen_ara = "Aarau"
# data management
ww_one = ww1 %>% 
  filter(ara_name==chosen_ara) %>% 
  # log
  mutate(logvl=log(vl)) %>%
  # duplicate indexes for INLA
  mutate(day1=day) %>% 
  mutate(below_loq=as.numeric(as.character(below_loq))) %>% 
  rownames_to_column()
# plot
g1 = ggplot(ww_one) +
  geom_point(aes(x=date,y=logvl),colour=cust_cols[1]) +
  labs(x="Date",y="Viral load (log2)")
g2 = ggplot(ww_one) +
  geom_histogram(aes(x=logvl),fill=cust_cols[1],colour="black",bins=30) +
  scale_y_continuous(expand=expansion(c(0,0.05))) +
  labs(x="Viral load (log2)",y="Count")
cowplot::plot_grid(g1,g2,labels=c("A","B"),rel_widths = c(2,1))

#'   
#' ## Model A1: linear regression
#' 
#' We use gamma regression as the dependent variable is continuous and strictly positive. The log transformation implies that all effects are multiplicative. We model the limit of quantification by allowing for a larger variability in viral load when the concentration estimate is below the LOQ (note that because of changing flow it may not always correspond to the lowest viral load values).
#' 
#' The first, basic model only includes a exponential trend over time:
#' 
#' $$
#' \log(V) = \alpha + \beta t
#' $$
#' 
ma1 = INLA::inla(vl ~ 1 + 
                   f(below_loq,model="iid") +
                   day,
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE,config=TRUE),
                 control.predictor = list(compute=TRUE,link=1))
summary(ma1)
summary_exp_vl(ma1,pars="day")
if(controls$compute_cv) kfoldcv_vl(ww_one,ma1)
ppp_vl(ww_one,ma1)

#' In the output we can see the call, the used time in seconds, the posterior estimates of the parameters and hyperparameters, and the WAIC. We show the exponential of the fixed parameters that are directly interpretable as multiplicative effects. We also look at the posterior predictive plot comparing the data with model predictions. 
#'   
#' ## Model A2: random walk
#' 
#' We model the variation in log viral load over time with a random walk, so that the difference between two successive observations follows a normal distribution.

ma2 = INLA::inla(vl ~ 1 +
                   f(below_loq,model="iid") +
                   f(day, model="rw1", scale.model=TRUE, constr=TRUE),
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE,config=TRUE),
                 control.predictor = list(compute=TRUE,link=1))
summary(ma2)
if(controls$compute_cv) kfoldcv_vl(ww_one,ma2)
ppp_vl(ww_one,ma2)

#' We see that the model fits much better, although visual inspection indicates over-fitting.
#' 
#' ## Model A3: random walk of order 2
#' 
#' We attempt to reduce over-fitting by using a random-walk of order 2, so that the difference depends on the last 2 observations.

ma3 = INLA::inla(vl ~ 1 +
                   f(below_loq,model="iid") +
                   f(day,model="rw2", scale.model=TRUE),
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE,config=TRUE),
                 control.predictor = list(compute=TRUE,link=1))
summary(ma3)
if(controls$compute_cv) kfoldcv_vl(ww_one,ma3)
ppp_vl(ww_one,ma3)

#' We still observe some over-fitting, even though there is a clear improvement. 
#' 
#' ## Model A4: random walk of order 2 + penalized complexity
#' 
#' We reduce over-fitting by adding priors that penalize complexity (Simpson et al, 2017).
#'  
ma4 = INLA::inla(vl ~ 1 +
                   f(below_loq,model="iid") +
                   f(day, model="rw2", scale.model=TRUE, constr=TRUE,
                     hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))),
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE,config=TRUE),
                 control.predictor = list(compute=TRUE,link=1))
summary(ma4)
if(controls$compute_cv) kfoldcv_vl(ww_one,ma4)
ppp_vl(ww_one,ma4)

#' The model now appropriately captures the dynamics of SARS-CoV-2 viral load, although still underestimates the variability especially during certain periods.
#'  
#' ## Model A5: additional covariates
#'  
#' There are other covariates available at the level of an ARA that may explain part of the residual variability. One is a change in the methods used to quantify SARS-CoV-2. A second is the day of the week, as behaviour may change for instance during weekends. For fixed effects we use weak normal priors with mean 0 and variance 5 on the log scale, corresponding to relative viral load values ranging between 0.01 and 80 (95% CrI).

ma5 = INLA::inla(vl ~ 1 +
                   f(below_loq,model="iid",constr=FALSE) +
                   f(day,model="rw2", scale.model=TRUE, constr=TRUE,
                     hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                   f(method,model="linear",mean.linear=0,prec.linear=.2) +
                   f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
                   f(hol,model="linear",mean.linear=0,prec.linear=.2),
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE,config=TRUE),
                 control.predictor = list(compute=TRUE,link=1))
summary(ma5)
summary_exp_vl(ma5,pars = "method|weekend|hol|test")
if(controls$compute_cv) kfoldcv_vl(ww_one,ma5)
ppp_vl(ww_one,ma5)


#' 
#' While the effect size is small, we observe an improvement in the fit. We also see the effect of method change, weekends and holidays.
#' 
#' We select this model as it gives the best compromise between accuracy (measured by RMSE), coverage and sharpness. We now apply this model to other ARAs.
#' 
#' 
#' ### Model A5 on Basel
#' 
#' We apply model A5 to the Basel ARA, and even though there is no change in method the model fits quite well.

chosen_ara = "Basel"
ww_one = ww1 %>% 
  filter(ara_name==chosen_ara) %>% 
  # log
  mutate(logvl=log(vl)) %>%
  # duplicate indexes for INLA
  mutate(day1=day) %>% 
  rownames_to_column()
ma5b = INLA::inla(vl ~ 1 +
                    f(below_loq, model="iid") +
                    f(day,model="rw2", scale.model=TRUE, constr=TRUE,
                      hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                    f(method,model="linear",mean.linear=0,prec.linear=.2) +
                    f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
                    f(hol,model="linear",mean.linear=0,prec.linear=.2),
                  data = ww_one,
                  family = "gamma",
                  control.compute = list(waic=TRUE,config=TRUE),
                  control.predictor = list(compute=TRUE,link=1))
summary(ma5b)
summary_exp_vl(ma5b,pars = "method|weekend|hol")
ppp_vl(ww_one,ma5b)


#' 
#' ### Model A5 on Region Bern
#' 
#' Contrary to Aarau, we find measurements of zero in Region Bern, which creates problems with logarithms. Instead of using zero-inflated models, we treat the zero values as very low (value = 1) and allow for additional variability, similarly to values below the limit of detection.
#' 

chosen_ara = "Region Bern"
# data management
ww_one = ww1 %>% 
  filter(ara_name==chosen_ara) %>% 
  # log
  mutate(logvl=log(vl)) %>%
  mutate(vl=if_else(vl==0, 1, vl)) %>%
  # mutate(below_loq=if_else(below_lod==1,below_lod,below_loq)) %>%
  # duplicate indexes for INLA
  mutate(day1=day) %>% 
  rownames_to_column()
ma5b = INLA::inla(vl ~ 1 +
                    f(below_loq,model="iid") +
                    f(below_lod,model="iid") +
                   f(day,model="rw2", scale.model=TRUE, constr=TRUE,
                     hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                    f(method,model="linear",mean.linear=0,prec.linear=.2) +
                    f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
                    f(hol,model="linear",mean.linear=0,prec.linear=.2),
                 data = ww_one,
                 family = "gamma",
                 control.compute = list(waic=TRUE,config=TRUE),
                 control.predictor = list(compute=TRUE,link=1))
summary(ma5b)
summary_exp_vl(ma5b,pars = "method|weekend|hol")
ppp_vl(ww_one,ma5b)



#' 
#' ### Model A5 on Saanen
#' 
#' Saanen is the place with the most number of measurements below the limit of detection.
#' 

chosen_ara = "Saanen"
# data management
ww_one = ww1 %>% 
  filter(ara_name==chosen_ara) %>% 
  # log
  filter(!is.na(vl)) %>% 
  mutate(logvl=log(vl)) %>%
  mutate(vl=if_else(vl==0, 1, vl)) %>%
  mutate(below_loq=if_else(below_lod==1,below_lod,below_loq)) %>% 
  # duplicate indexes for INLA
  mutate(day1=day) %>% 
  rownames_to_column()
ma5b = INLA::inla(vl ~ 1 +
                    f(below_loq,model="iid") +
                    f(below_lod,model="iid") +
                    f(day,model="rw2", scale.model=TRUE, constr=TRUE,
                      hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                    f(method,model="linear",mean.linear=0,prec.linear=.2) +
                    f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
                    f(hol,model="linear",mean.linear=0,prec.linear=.2),
                  data = ww_one,
                  family = "gamma",
                  control.compute = list(waic=TRUE,config=TRUE),
                  control.predictor = list(compute=TRUE,link=1))
summary(ma5b)
summary_exp_vl(ma5b,pars = "method|weekend|hol")
ppp_vl(ww_one,ma5b)

#' 
#' ### Model A5 on Laupen
#' 
#' In Laupen there are many measurements below the limit of quantification, the model is still fine.
#' 

chosen_ara = "Laupen"
# data management
ww_one = ww1 %>% 
  filter(ara_name==chosen_ara) %>% 
  # log
  mutate(logvl=log(vl)) %>%
  mutate(vl=if_else(vl==0, 1, vl)) %>%
  # duplicate indexes for INLA
  mutate(day1=day) %>% 
  rownames_to_column()
ma5b = INLA::inla(vl ~ 1 +
                    f(below_loq,model="iid") +
                    f(below_lod,model="iid") +
                    f(day,model="rw2", scale.model=TRUE, constr=TRUE,
                      hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                    f(method,model="linear",mean.linear=0,prec.linear=.2) +
                    f(weekend,model="linear",mean.linear=0,prec.linear=.2) +
                    f(hol,model="linear",mean.linear=0,prec.linear=.2),
                  data = ww_one,
                  family = "gamma",
                  control.compute = list(waic=TRUE,config=TRUE),
                  control.predictor = list(compute=TRUE,link=1))
summary(ma5b)
summary_exp_vl(ma5b,pars = "method|weekend|hol")
ppp_vl(ww_one,ma5b)
