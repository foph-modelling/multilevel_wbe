

true_shapes = fread('data/true_shapes.csv')
true_shapes = true_shapes[period %in% c(1,2),]
true_shapes[, mean_disc := infotheo::discretize(mean, )]



entropies = data.table(ara_id=unique(ww_all$ara_id), 
                      entropy= 
                        sapply(unique(ww_all$ara_id), 
                               function(x){infotheo::entropy(true_shapes[ara_id == x]$mean_disc)})
                      )

mi_old=100
for(i in 1:1000){
  aras = sample(unique(ww_all$ara_id), 10, replace = F)
  mi = infotheo::multiinformation(tail(dcast(true_shapes[ara_id %in% aras], day ~ ara_id, value.var = 'mean_disc'), -6)[,-c('day')])
  nmi = mi / sum(entropies[ara_id %in% aras,]$entropy)
  if(mi < mi_old){
    mi_old = nmi
    aras_select = aras
  }
  
}



  
  
  func1 = function(n){
  
  
  func_nmi = function(x){
    aras = sample(unique(ww_all$ara_id), n, replace = F)
    mi = infotheo::multiinformation(
      tail(
        dcast(
          true_shapes[ara_id %in% aras], 
          day ~ ara_id, 
          value.var = 'mean_disc'), 
        -6)[,-c('day')])
    nmi = mi / sum(entropies[ara_id %in% aras,]$entropy)
    nmi}
  
  mis = sapply(1:1000, func)
  mean(mis)
  
  }
  
  mean_mis = sapply(seq(20,1,-1),func1)

  
func_nmi = function(x){
    aras = sample(unique(ww_all$ara_id), n, replace = F)
    mi = infotheo::multiinformation(
      tail(
        dcast(
          true_shapes[ara_id %in% aras], 
          day ~ ara_id, 
          value.var = 'mean_disc'), 
        -6)[,-c('day')])
    nmi = mi / sum(entropies[ara_id %in% aras,]$entropy)
    nmi}
  
func_ent = function(x){
  aras = sample(unique(ww_all$ara_id), n, replace = F)
  mi = infotheo::entropy(
    tail(
      dcast(
        true_shapes[ara_id %in% aras], 
        day ~ ara_id, 
        value.var = 'mean_disc'), 
      -6)[,-c('day')]) / infotheo::entropy(
        tail(
          dcast(
            true_shapes, 
            day ~ ara_id, 
            value.var = 'mean_disc'), 
          -6)[,-c('day')])
  tui = mi
  tui}


tis_dt = data.table()
mis_dt = data.table()
for(n in seq(20,2,-1)){
  
    mis = sapply(1:1000, func_nmi)
    
    mis_dt = rbind(
      mis_dt, 
      data.table(n=n, mean = mean(mis), upper = quantile(mis, 0.975), lower=quantile(mis, 0.025))
    )
    
    tis = sapply(1:1000, func_ent)
    tis_dt = rbind(
      tis_dt, 
      data.table(n=n, mean = mean(tis), upper = quantile(tis, 0.975), lower=quantile(tis, 0.025))
    )
    
}

tis_dt[, type := 'Proportion of total information conserved']
mis_dt[, type := 'Normalised mutual information']

information_dt = rbind(tis_dt, mis_dt)

information_dt %>% ggplot() + 
  geom_pointrange(aes(x=n, y=mean, ymin=lower, ymax=upper)) + 
  facet_wrap(~type, ncol=1, scales='free') + 
  theme_minimal()


n = 4

nmi_tui = data.table()
for(i in 1:1000){
  aras = sample(unique(ww_all$ara_id), n, replace = F)
  mi = infotheo::multiinformation(
    tail(
      dcast(
        true_shapes[ara_id %in% aras], 
        day ~ ara_id, 
        value.var = 'mean_disc'), 
      -6)[,-c('day')])
  nmi = mi / sum(entropies[ara_id %in% aras,]$entropy)
  nmi

  tui = infotheo::entropy(
    tail(
      dcast(
        true_shapes[ara_id %in% aras], 
        day ~ ara_id, 
        value.var = 'mean_disc'), 
      -6)[,-c('day')]) / infotheo::entropy(
        tail(
          dcast(
            true_shapes, 
            day ~ ara_id, 
            value.var = 'mean_disc'), 
          -6)[,-c('day')])
  
  nmi_tui = rbind(nmi_tui, data.table(nmi = nmi, tui = tui, mi = mi))}

nmi_tui %>% ggplot() + geom_point(aes(x=nmi, y=tui))


mi_old=0
for(i in 1:1000){
  aras = sample(unique(ww_all$ara_id), 6, replace = F)
  mi = infotheo::multiinformation(tail(dcast(true_shapes[ara_id %in% aras], day ~ ara_id, value.var = 'mean_disc'), -6)[,-c('day')])
  nmi = mi / sum(entropies[ara_id %in% aras,]$entropy)
  
  tui = infotheo::entropy(
    tail(
      dcast(
        true_shapes[ara_id %in% aras], 
        day ~ ara_id, 
        value.var = 'mean_disc'), 
      -6)[,-c('day')]) / infotheo::entropy(
        tail(
          dcast(
            true_shapes, 
            day ~ ara_id, 
            value.var = 'mean_disc'), 
          -6)[,-c('day')])
  
  
  if(tui > mi_old){
    mi_old = tui
    aras_select = aras
  }
  
}


mi_old
aras_select

catchments %>% 
  filter(ara_id %in% aras_select) %>% 
  ggplot() + geom_sf(data=shapes$canton_shp) + 
  geom_sf(aes(fill=ara_id))
