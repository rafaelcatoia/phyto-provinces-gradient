#list_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_geo_abiotics_dists'))


#function evaluate_cluster
## this function runs inside eval_clustering_OOS function.
## it receives a distance Matrix that will be evaluated
## and a data frame containg the cluster membership of each samples.

## Objects necessary for debugging it::
#distMatrixEval = list_geo_abiotics_dists$geo_Dist_mirrored
#list_geo_abiotics_dists$geo_Dist %>% plot_distance_matrix
#list_geo_abiotics_dists$geo_Dist_mirrored %>% plot_distance_matrix
#list_geo_abiotics_dists$abioticDist %>% plot_distance_matrix
#list_geo_abiotics_dists$aitDist %>% plot_distance_matrix
evaluate_cluster <- function(distMatrixEval,df_clusters,Kmax=Kmax){
  out_df=data.frame(
    ncluster = 0,
    ## Within ---------
    within_sum = 0,
    within_sum_sq = 0 ,
    max_max_dist_within = 0,
    avg_max_dist_within = 0,
    avg_avg_within_pairwise_dist = 0,
    ## between ---------
    avg_dist_between_medoids = 0,
    avg_dist_within_dist2medoid = 0,
    ## Ratio -----------
    avg_within_between_dist2medoid = 0
  )
  
  for(itt in 2:Kmax){
    dfX = data.frame(X=df_clusters[,itt])
    numberOfClusters_itt <- max(dfX$X)
    vet_number_obs_per_clusters <- dfX %>% group_by(X) %>% summarise(n()) %>% pull()
    idx_list <- list()
    medoid_idx = rep(NA,length = numberOfClusters_itt)
    medoid_idx_within = rep(NA,length = numberOfClusters_itt)
    distMatrix_subseted <- list()
    
    for(iii in 1:numberOfClusters_itt){
      # getting the index of each observation for each clusters
      idx_clust_iii = which(dfX$X==iii)
      
      # getting the index of each observation for each clusters --- storing it
      idx_list[[iii]] <- idx_clust_iii
      
      # subsetting the distance matrix
      distMatrix_subseted[[iii]] = distMatrixEval[idx_clust_iii,idx_clust_iii]
      
      #getting the medoid id for each cluster 
      #first we store the subsetted idx for the medoid
      if(length(idx_clust_iii)<2){
        medoid_idx_within[iii] = 1
      }else{
        medoid_idx_within[iii] = cluster::pam(x = as.dist(distMatrix_subseted[[iii]]),k = 1)$id.med
      }
      #than we store the idx in the not subsetted version.
      medoid_idx[iii] = idx_clust_iii[medoid_idx_within[iii]]
      # cat('\n Done for cluster ', i,'---------------------- \n')
    }
    
    ## Sum within pairwise dist
    sum_within_pairwise_dist         = lapply(distMatrix_subseted,function(x) {sum(x)/2}) %>% unlist() %>% sum()
    sum_within_pairwise_dist_squared = lapply(distMatrix_subseted,function(x) {((x^2 %>% sum())/2)}) %>% unlist() %>% sum()
    avg_within_pairwise_dist = lapply(distMatrix_subseted,function(x) { if(!is.null(dim(x))){(sum(x)/2)/sum(lower.tri(x)) }else{0}}) %>% unlist %>% mean
    ## Metrics within -------------------------------------------------------------------------
    sum_within_dist2medoid <- rep(NA,length = numberOfClusters_itt)
    avg_within_dist2medoid <- rep(NA,length = numberOfClusters_itt)
    max_within_dist2medoid <- rep(NA,length = numberOfClusters_itt)
    max_within_dist <- rep(NA,length = numberOfClusters_itt)
    min_within_dist <- rep(NA,length = numberOfClusters_itt)
    
    for(i in 1:numberOfClusters_itt){
      if(is.null(dim(distMatrix_subseted[[i]]))){
        #value if the cluster is composed by only one observation
        sum_within_dist2medoid[i] = 0
        avg_within_dist2medoid[i] = 0
        max_within_dist2medoid[i] = 0
        max_within_dist[i]=0
        min_within_dist[i]=0
      }else{
        sum_within_dist2medoid[i] = sum(distMatrix_subseted[[i]][medoid_idx_within[i],])
        avg_within_dist2medoid[i] = sum_within_dist2medoid[i]/(nrow(distMatrix_subseted[[i]])-1)
        max_within_dist2medoid[i] = max(distMatrix_subseted[[i]][medoid_idx_within[i],])
        max_within_dist[i]=max(distMatrix_subseted[[i]])
        min_within_dist[i]=min(distMatrix_subseted[[i]][upper.tri(distMatrix_subseted[[i]])])
      }
    }
    
    ## Metrics between -------------------------------------------------------------------------
    combinations_between_clusters = choose(numberOfClusters_itt,2)
    sum_dist_between_medoid <- rep(NA,length = numberOfClusters_itt)
    avg_dist_between_medoid <- rep(NA,length = numberOfClusters_itt)
    max_dist_between_medoid <- rep(NA,length = numberOfClusters_itt)
    min_dist_between_medoid <- rep(NA,length = numberOfClusters_itt)
    max_dist_between_dist <- rep(NA,length = numberOfClusters_itt)
    min_dist_between_dist <- rep(NA,length = numberOfClusters_itt)
    
    ### Calculating  ---------------------------------------------------------------------------
    dist_matrix_medoids = distMatrixEval[medoid_idx,medoid_idx]
    
    ## sum of distances between medoids
    sum_dist_between_medoid = sum(dist_matrix_medoids)/2
    sum_dist_between_medoid_squared =  (sum(dist_matrix_medoids)/2)^2
    ## avg of distances between medoids
    avg_dist_between_medoid = sum_dist_between_medoid/combinations_between_clusters
    avg_dist_between_medoid_squared = sum_dist_between_medoid/combinations_between_clusters
    ## max and min between medoids
    max_dist_between_medoid = max(dist_matrix_medoids)
    min_dist_between_medoid = min(dist_matrix_medoids[upper.tri(dist_matrix_medoids)])
    
    ## max between clusters and min between clusters
    min_dist_between_clusters = matrix(0,nrow = numberOfClusters_itt,ncol=numberOfClusters_itt)
    for(i in 1:numberOfClusters_itt){
      for(j in 1:numberOfClusters_itt){
        if(i>j){
          dist_aux = distMatrixEval[idx_list[[i]],idx_list[[j]]]
          min_dist_between_clusters[i,j] <- min(dist_aux)
        }
      }
    }
    
    min_dist_between_clusters=min_dist_between_clusters[lower.tri(min_dist_between_clusters)]
    avg_min_dist_between_clusters = mean(min_dist_between_clusters)
    sum_min_dist_between_clusters = sum(min_dist_between_clusters)
    
    ## arranging output
    out <- data.frame(
      ncluster = numberOfClusters_itt,
      ## Within ---------
      within_sum = sum_within_pairwise_dist            ,
      within_sum_sq = sum_within_pairwise_dist_squared ,
      max_max_dist_within = max(max_within_dist)       ,
      avg_max_dist_within = mean(max_within_dist)      ,
      avg_avg_within_pairwise_dist = avg_within_pairwise_dist,
      ## between ---------
      avg_dist_between_medoids = avg_dist_between_medoid,
      avg_dist_within_dist2medoid = mean(avg_within_dist2medoid),
      ## Ratio -----------
      avg_within_between_dist2medoid = (mean(avg_within_dist2medoid) / avg_dist_between_medoid)
    )
    out_df = bind_rows(out_df,out)
  }
  out_df = out_df %>% filter(ncluster>1)
  return(out_df)
}







