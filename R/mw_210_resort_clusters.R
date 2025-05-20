library(dplyr)

# Sample function to calculate Jaccard similarity between two sets
jaccard_similarity <- function(set1, set2) {
  length(intersect(set1, set2))
}



tt_cluster_select = tt_cluster_all %>% filter((period==1 & nclust==7) | (period==2 & nclust==5) | (period==3 & nclust==4) | (period==4 & nclust == 2) | (period == 5 & nclust==3) ) 

find_cluster_swaps = function(tt_cluster_select){
  swap_from_all = list()
  swap_to_all = list()
  for (peri in 2:5){
    
    curr_part = peri
    prev_part = 1
    clusts_in_curr_part = sort(unique((tt_cluster_select %>% filter(period == curr_part))$cluster))
    clusts_in_prev_part = sort(unique((tt_cluster_select %>% filter(period == prev_part))$cluster))
    nrow = length(clusts_in_curr_part)
    ncol = length(clusts_in_prev_part)
    sim_mat = matrix(ncol = ncol, nrow=nrow)
    for (i in 1:length(clusts_in_curr_part)){
      for (j in 1:length(clusts_in_prev_part)){
        
        print(i,j)
        
        sim_mat[i,j] = jaccard_similarity( unique((tt_cluster_select %>% filter(period == curr_part, cluster == i))$ara_name),  unique((tt_cluster_select %>% filter(period == prev_part, cluster == j))$ara_name))
        
      }
    }
    
    swap_from = c()
    swap_to = c()
    for(clust in 1:length(clusts_in_curr_part)){
      
      max_val =max(sim_mat)
      
      if(max_val>0){
        ind = which(sim_mat== max_val)[1]
        
        i = ind%%nrow
        j = ind%/%nrow
        
        if(i==0){
          i=nrow
        }else{
          j=j+1
        }
        swap_to = c(swap_to, j)
        swap_from = c(swap_from, i)
        
        print(ind%/%nrow )
        print(sim_mat[,j])
        print(sim_mat[i,])
        sim_mat[,j] = 0 
        sim_mat[i,] = 0 
      }else{
        i = setdiff(1:nrow, swap_from)[1]
        j = setdiff(1:ncol, swap_to)[1]
        swap_to = c(swap_to, j)
        swap_from = c(swap_from, i)
      }
      
        
      
      
    }
    swap_from_all[[peri]] = swap_from
    swap_to_all[[peri]] = swap_to
  }
  
  return(list(swap_from_all, swap_to_all))
}
