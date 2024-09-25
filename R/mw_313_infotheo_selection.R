true_shapes = fread('data/true_shapes.csv')
true_shapes = true_shapes[period %in% c(1,2),]
true_shapes[, mean_disc := infotheo::discretize(mean, nbins = 10)]

true_shapes %>% ggplot() + geom_point(aes(x=date, y=mean_disc)) + facet_wrap(~ara_id)




#using TS
all_nmis = data.table()

for(nsamp in seq(20, 2, -2)){
  nmis = data.table()
  for(i in 1:1000){
    aras = sample(unique(ww_all$ara_id), nsamp, replace = F)
    Y = tail(dcast(true_shapes[ara_id %in% aras,], day ~ ara_id, value.var = 'mean_disc'), -6)[,-c('day')]
    X = tail(dcast(true_shapes, day ~ ara_id, value.var = 'mean_disc'), -6)[,-c('day')]
    
    
    # Compute individual and joint entropies
    H_X <- entropy(X)
    H_Y <- entropy(Y)
    H_XY <- entropy(cbind(X, Y))
    H_XconY <- condentropy(X, Y)
    
    # Compute mutual information
    MI_XY <- mutinformation(X, Y)
    
    # Max-entropy normalization
    NMI_max_entropy <- MI_XY / max(H_X, H_Y)
    
    # Joint-entropy normalization
    NMI_joint_entropy <- MI_XY / H_XY
    
    
    # Joint-entropy normalization
    NMI_dij_entropy <- MI_XY / (H_X/100 + H_Y/nsamp)
    
    prop_information = 1 - H_XconY/H_X
    
    TC_X = multiinformation(X) / sum(sapply(1:ncol(X), 
                                            function(i){infotheo::entropy(X[,..i])}))#nsamp#entropy(X)
    TC_Y = multiinformation(Y) / sum(sapply(1:ncol(Y), 
                                            function(i){infotheo::entropy(X[,..i])}))#nsamp#entropy(Y) 
    
    ratio_of_total_correlation = TC_X / TC_Y
    
    
    # Print the results
    
    nmis = rbind(nmis, data.table(NMI_joint_entropy, prop_information, NMI_max_entropy, ratio_of_total_correlation, nsamp))
  }
  all_nmis = rbind(all_nmis, nmis)
  
}

Y = selected_sites

(multiinformation(X) / sum(sapply(1:ncol(X), function(i){infotheo::entropy(X[,..i])})) /(multiinformation(Y) / sum(sapply(1:ncol(Y), function(i){infotheo::entropy(Y[,..i])}))))

all_nmis %>% ggplot() + 
  geom_point(aes(x=NMI_joint_entropy, y=1./ratio_of_total_correlation, color=as.character(nsamp)), alpha=0.4) + 
  geom_point(x=entropy(selected_sites)/entropy(X), y= 1 / ( (multiinformation(X) / sum(sapply(1:ncol(X), function(i){infotheo::entropy(X[,..i])})) /(multiinformation(Y) / sum(sapply(1:ncol(Y), function(i){infotheo::entropy(Y[,..i])}))))) )+ 
  xlab('Proportion of information retained') + 
  ylab('Ratio of total correlation (log)') + 
  #scale_y_continuous(trans='log10')+
  #scale_x_continuous(trans='log10')+
  theme_bw() + 
  theme(
    legend.position = 'none'
  )


# Greedy algorithm for subset selection
subset_size <- 6
selected_indices <- c()
remaining_indices <- sample(1:ncol(X), ncol(X), replace = F)

for (i in 1:subset_size) {
  best_entropy_gain <- -Inf
  best_index <- NULL
  
  # Try adding each remaining time-series and compute joint entropy
  for (idx in remaining_indices) {
    print(idx)
    
    new_indicies = c(selected_indices, idx)
    current_subset <- X[,..new_indicies]
    current_entropy <- entropy(current_subset)
    
    if (current_entropy > best_entropy_gain) {
      best_entropy_gain <- current_entropy
      best_index <- idx
    }
  }
  
  # Add the best time-series to the selected subset
  selected_indices <- c(selected_indices, best_index)
  remaining_indices <- setdiff(remaining_indices, best_index)
}

# The selected_indices vector now contains the 6 selected time-series indices
print(selected_indices)


selected_sites = X[,..selected_indices]

entropy(selected_sites)/entropy(X)

selected_aras = colnames(selected_sites)

catchments %>% filter(ara_id %in% selected_aras) %>%  
    ggplot() + 
      geom_sf(data=shapes$canton_shp, fill=NA) + 
      geom_sf(data=shapes$see_shp, fill='lightblue', color=NA) + 
      geom_sf( fill= 'red') +
      theme_minimal()





















###### - - - - - - - - - - - -  - - - - - 


#using TS
all_nmis = data.table()

for(nsamp in seq(20, 2, -2)){
  nmis = c()
  for(i in 1:1000){
    aras = sample(unique(ww_all$ara_id), nsamp, replace = F)
    Y = tail(dcast(true_shapes[ara_id %in% aras,], day ~ ara_id, value.var = 'mean_disc'), -6)[,-c('day')]
    X = tail(dcast(true_shapes, day ~ ara_id, value.var = 'mean_disc'), -6)[,-c('day')]
    
    
    # Compute individual and joint entropies
    H_X <- entropy(X)
    H_Y <- entropy(Y)
    H_XY <- entropy(cbind(X, Y))
    
    # Compute mutual information
    MI_XY <- mutinformation(X, Y)
    
    # Max-entropy normalization
    NMI_max_entropy <- MI_XY / max(H_X, H_Y)
    
    # Joint-entropy normalization
    NMI_joint_entropy <- MI_XY / H_XY
    
    # Print the results
    
    nmis = c(nmis, NMI_max_entropy)
  }
  all_nmis = cbind(all_nmis, nmis)
  
}

colnames(all_nmis) = paste0(rep('samples_', 10), seq(20, 2, -2))


all_nmis_long = melt(all_nmis, measure.vars = colnames(all_nmis), variable.name = 'samples', value.name = 'NMI') 

all_nmis_long %>% ggplot() + 
  geom_density(aes(x=NMI, color=samples, group=samples))

all_nmis_long_both[samples != 'samples_2',] %>% ggplot() + 
  geom_density(aes(x=NMI, color=method, group=method)) + 
  facet_wrap(~samples)


# Using residuals 
all_nmis = data.table()

for(nsamp in seq(10, 2, -2)){
  nmis = c()
  for(i in 1:1000){
    aras = sample(unique(ww_all$ara_id), nsamp, replace = F)
    Y = tail(dcast(true_shapes[ara_id %in% aras,], day ~ ara_id, value.var = 'mean'), -6)[,-c('day')]
    X = tail(dcast(true_shapes, day ~ ara_id, value.var = 'mean'), -6)[,-c('day')]
    
    var_X <- vars::VAR(X, p = 1)  # p is the lag order
    var_Y <- vars::VAR(Y, p = 1)
    
    # Extract residuals after fitting VAR models
    residuals_X <- residuals(var_X)
    residuals_Y <- residuals(var_Y)
    
    
    # Discretize the residuals to prepare for MI calculation
    residuals_X_disc <- discretize(residuals_X, disc = "equalwidth", nbins = 100)
    residuals_Y_disc <- discretize(residuals_Y, disc = "equalwidth", nbins = 100)
    
    # Compute individual and joint entropies
    H_X <- entropy(residuals_X_disc)
    H_Y <- entropy(residuals_Y_disc)
    H_XY <- entropy(cbind(residuals_X_disc, residuals_Y_disc))
    
    # Compute mutual information
    MI_XY <- mutinformation(residuals_X_disc, residuals_Y_disc)
    
    # Max-entropy normalization
    NMI_max_entropy <- MI_XY / max(H_X, H_Y)
    
    # Joint-entropy normalization
    NMI_joint_entropy <- MI_XY / H_XY
    
    
    # Joint-entropy normalization
    NMI_dij_entropy <- MI_XY / (H_X/100 + H_Y/nsamp)
    
    # Print the results
    
    nmis = c(nmis, NMI_max_entropy)
  }
  all_nmis = cbind(all_nmis, nmis)
  
}

colnames(all_nmis) = paste0(rep('samples_', 5), seq(10, 2, -2))
all_nmis_long = melt(all_nmis, measure.vars = colnames(all_nmis), variable.name = 'samples', value.name = 'NMI') 

all_nmis_long %>% ggplot() + 
  geom_histogram(aes(x=NMI, fill=samples, group=samples))

entropies = data.table(ara_id=unique(true_shapes$ara_id), 
                       entropy= 
                         sapply(unique(true_shapes$ara_id), 
                                function(x){infotheo::entropy(true_shapes[ara_id == x]$mean_disc)})
)


colnames(all_nmis) = paste0(rep('samples_', 10), seq(20, 2, -2))


all_nmis_long_ce = melt(all_nmis, measure.vars = colnames(all_nmis), variable.name = 'samples', value.name = 'NMI') 

all_nmis %>% ggplot() + 
  geom_density(aes(y=prop_information, color=as.character(nsamp), group=nsamp))


all_nmis %>% ggplot() + 
  geom_violin(aes(x=nsamp, y=prop_information, group=nsamp, color=as.character(nsamp)))+ 
  ylab('proportion of information retained') + 
  xlab('number of treatment plants sampled') + 
  theme(
    legend.position = 'none'
  )



#using TS
all_nmis = data.table()

for(nsamp in seq(20, 2, -2)){
  nmis = c()
  for(i in 1:1000){
    aras = sample(unique(ww_all$ara_id), nsamp, replace = F)
    Y = tail(dcast(true_shapes[ara_id %in% aras,], day ~ ara_id, value.var = 'mean_disc'), -6)[,-c('day')]
    X = tail(dcast(true_shapes, day ~ ara_id, value.var = 'mean_disc'), -6)[,-c('day')]
    
    
    # Compute individual and joint entropies
    H_X <- entropy(X)
    H_Y <- entropy(Y)
    H_XY <- entropy(cbind(X, Y))
    H_XconY <- condentropy(X, Y)
    
    prop_information = sum(sapply(1:nsamp, FUN = function(x){mutinformation(Y[,..x], X)}))/
    
    
    # Print the results
    
    nmis = c(nmis, prop_information)
  }
  all_nmis = cbind(all_nmis, nmis)
  
}


colnames(all_nmis) = paste0(rep('samples_', 10), seq(20, 2, -2))


all_nmis_long_si = melt(all_nmis, measure.vars = colnames(all_nmis), variable.name = 'samples', value.name = 'NMI') 

all_nmis_long_si %>% ggplot() + 
  geom_density(aes(x=NMI, color=samples, group=samples))


all_nmis_long_si %>% ggplot() + 
  geom_violin(aes(x=samples, y=nmi_norm, group=samples, color=samples))





# using a target variable PCA based 
X = tail(dcast(true_shapes, day ~ ara_id, value.var = 'mean_disc'), -6)[,-c('day')]

pca <- prcomp(X)

# Use the first principal component as the target variable (this captures the most variance)
Z <- pca$x[,1]  # First principal component


Z <- discretize(as.data.frame(Z), nbins=10)


all_nmis = data.table()

for(nsamp in seq(20, 2, -2)){
  nmis = data.table()
  for(i in 1:1000){

    aras = sample(unique(ww_all$ara_id), nsamp, replace = F)
    Y = tail(dcast(true_shapes[ara_id %in% aras,], day ~ ara_id, value.var = 'mean_disc'), -6)[,-c('day')]
    X = tail(dcast(true_shapes, day ~ ara_id, value.var = 'mean_disc'), -6)[,-c('day')]
    Y_prime = tail(dcast(true_shapes[!(ara_id %in% aras),], day ~ ara_id, value.var = 'mean_disc'), -6)[,-c('day')]
    
    # 1. Calculate Mutual Information between X10 and Y (Total Information in X10 about Y)
    MI_Y_Z <- mutinformation(Y, Z)
    
    MI_Y_prime_Z <- mutinformation(Y_prime, Z)
    
    
    # 2. Calculate Mutual Information between X100 and Y (Total Information in X100 about Y)
    MI_X_Z <- mutinformation(X, Z)
    
    # 3. Calculate Redundant Information:
    # The redundancy can be estimated as the minimum mutual information between X10 and other variables in X100 with respect to Y
    # For simplicity, we'll take the minimum MI between each variable in X100 and Y as a proxy
    #MI_X_vars_Z <- apply(X, 2, function(var) mutinformation(as.data.frame(var), Z))
    Redundant_info <- min(MI_Y_Z, MI_X_Z)
    
    # 4. Calculate Unique Information:
    # Unique Information is the difference between the total mutual information in X10 and the redundant information
    Unique_info_Y <- MI_Y_Z - Redundant_info
    
    Unique_info_Y_prime <- MI_Y_prime_Z - Redundant_info
    
    Unique_info_X <- MI_X_Z - Redundant_info
    
    # 5. Calculate Synergistic Information:
    # Synergy can be seen as information not captured by the individual parts but only by combining X10 and X100
    # Synergistic information can be approximated as:
    Synergistic_info <- MI_Y_Z - MI_Y_prime_Z - Unique_info_X
    
    # Print the results
    #cat("Mutual Information (X10, Y):", MI_Y_Z, "\n")
    #cat("Mutual Information (X100, Y):", MI_X_Z, "\n")
    #cat("Redundant Information:", Redundant_info, "\n")
    #cat("Unique Information:", Unique_info_Y_prime, "\n")
    #cat("Synergistic Information:", Synergistic_info, "\n")

    nmis = rbind(nmis, data.table(Redundant_info, Unique_info_Y, Unique_info_Y_prime, Unique_info_X, nsamp))
  }
  all_nmis = rbind(all_nmis, nmis)
  
}

condinformation(Y_prime, Z, Y)

all_nmis %>% ggplot() + 
  geom_violin(aes(x=nsamp, y=Unique_info_Y/max(Unique_info_Y), group=nsamp, color=nsamp))
