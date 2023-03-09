# multilevel_wbe
Bayesian multilevel models for wastewater-based epidemiology

Model we could fit to the data using INLA in R (https://www.r-inla.org/):

i = ARA
j = testing lab
k = time period

log(conc_{ijk}(t_k)) = alpha_i + beta_{ik} t_k + gamma_j + delta v + noise

Where alpha_i is the effect of the ARA, beta_{ik} is the slope on time period t_k, gamma_j is the lab specific effect and delta the variant effect. 

# Open questions

* Should we take the log of the log(concentration) data? 
* How do we want to include the variant effect? Using the relative fraction of each variant? Would we have enough data to include this, i.e. for each location and time point an estimate? Or can we only do this is segregated data. 
* Could the testing lab also influence the slope? Or even the slope rather than the intercept?
* Same for the variant, this effect maybe needs to include a time effect. 
