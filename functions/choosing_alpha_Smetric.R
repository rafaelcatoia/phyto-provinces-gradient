############# -------------------------------------------
## script to create a function to evaluate the clustering
############# -------------------------------------------
 
# df_geo_abiotics = readRDS(file = paste0(savingdir,'/','df_geo_abiotics'))
# df_geo_abiotics = df_geo_abiotics %>% arrange(SampleID)
# list_normalized_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_normalized_geo_abiotics_dists'))
# list_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_geo_abiotics_dists'))
# list_AitDist = readRDS(file = paste0(savingdir,'/','list_AitDist'))
# D = list_AitDist[[1]]

# list_geo_abiotics_dists$geo_Dist_mirrored %>% plot_distance_matrix()
# list_geo_abiotics_dists$abioticDist %>% plot_distance_matrix()
# 

choosing_alpha_Smetric <- function(
  D,
  list_normalized_dist=list_normalized_geo_abiotics_dists, # those dists can be used to create a convex mixture of distance matrices
  list_dist = list_geo_abiotics_dists, # those will be use to perform cluser evaluation
  latMirrored=T,grid_alpha=seq(0,1,0.1),
  K=10){
  
  
  nrows = dim(D$INS)[1]
  
  mat_hclust <- matrix(NA,nrow=nrows,ncol=length(grid_alpha))
  #mat_pam <- matrix(NA,nrow=nrows,ncol=length(vet_alpha))
  list_df_eval = list()
  for( i in 1:length(grid_alpha)){
    D_mix = (1-grid_alpha[i]) * D$INS + grid_alpha[i]*list_normalized_dist$geoDist
    mat_hclust[,i] <- cutree(hclust(d=as.dist(D_mix),method='ward.D'),k = K)
    #mat_pam[,i] <- cluster::pam(x=as.dist(D_mix),k = K,cluster.only = T)
    
  
    list_df_eval[[i]]<- bind_rows(
      clustEvaluation(distMatrixEval = D$OOS,vetClusterMembership = mat_hclust[,i]) %>%
        mutate(method='hclust',alpha=grid_alpha[i], DistMetric='Biotic') ,
      
      clustEvaluation(distMatrixEval = list_dist$abioticDist,vetClusterMembership = mat_hclust[,i]) %>%
        mutate(method='hclust',alpha=grid_alpha[i], DistMetric='Abiotic'), 
      
      clustEvaluation(distMatrixEval = list_dist$geoDist, vetClusterMembership = mat_hclust[,i]) %>%
        mutate(method='hclust',alpha=grid_alpha[i], DistMetric='Geo')
      #compute_medoid_distances(distance_matrix = list_dist$geo_Dist,
      #                         group_vector = mat_pam[,i]) %>% mutate(method='pam',alpha=vet_alpha[i], DistMetric='Geo'), 
      #compute_medoid_distances(distance_matrix = list_dist$geo_Dist_mirrored,
      #                         group_vector = mat_pam[,i]) %>% mutate(method='pam',alpha=vet_alpha[i], DistMetric='Geo_Mirrored'), 
      #compute_medoid_distances(distance_matrix = list_dist$abioticDist,
      #                         group_vector = mat_pam[,i]) %>% mutate(method='pam',alpha=vet_alpha[i], DistMetric='Abiotics'), 
      #compute_medoid_distances(distance_matrix = list_dist$aitDist,
      #                         group_vector = mat_pam[,i]) %>% mutate(method='pam',alpha=vet_alpha[i], DistMetric='Aitchison')
    )
  }
  df_out <- data.table::rbindlist(list_df_eval)
  return(df_out) 
}
