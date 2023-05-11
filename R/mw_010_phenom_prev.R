#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: phenomenological model linking concentration to prevalence from Morvan et al (Nat Com 2022)
# creation: jriou & mwagner
# init date: 2023-04-17
#:::::::::::::::::::::::::::::

mw_010_phenom_prev = function(dat) {
  
  # Qp: wastewater generated per person per day (L/day)
  Qp = 400
  # S: mean shedding rate (gc/mL)
  S = 1.9e6
  # V: mean volume of stool per person per day (mL/day)
  V = 128*1.06
  
  res = dat %>% 
    dplyr::mutate(prev = conc * Qp / (S*V))
  res %>% 
    filter(ara_name=="Aarau") %>% 
    ggplot() +
    # geom_point(aes(x=date,y=conc)) +
    geom_point(aes(x=date,y=prev),col='red')
}