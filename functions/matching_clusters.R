##### Matching Function 

## objects for debugging 
#DistMatrix = list_dist_toMatch$geoDist
#trueLabel = trueClusterMembership[,i]
#proposedLabel = mat_cluster_membership[,i]

matching_clusters <- function(DistMatrix,trueLabel,proposedLabel){

  Kmatch = max(trueLabel)
  matrixMatching_mean <- matrix(NA,nrow=Kmatch,ncol=Kmatch)
  #matrixMatching_max <- matrix(NA,nrow=Kmatch,ncol=Kmatch)
  
  for(ii in 1:Kmatch){ 
    for(jj in 1:Kmatch){
      idx_trueLabel = which(trueLabel==ii)
      idx_replicateLabel = which(proposedLabel==jj)
      subsetedDistMatrix <- DistMatrix[idx_trueLabel,idx_replicateLabel]
      matrixMatching_mean[[ii,jj]]<- mean(subsetedDistMatrix)
      #matrixMatching_max[[ii,jj]] <- batata(subsetedDistMatrix)
      #cat(ii,jj,'\n')
    }
  }
  
  ord_mean = clue::solve_LSAP(matrixMatching_mean) %>% as.numeric()
  #ord_max = clue::solve_LSAP(matrixMatching_max) %>% as.numeric()
  
  matchedMean = relabel_matching(x = proposedLabel,ord_mean)
  #matchedMax = relabel_matching(x = proposedLabel,ord_max)
  
  return(matchedMean)
  #list(
    #matchedMean = matchedMean#,
    #matchedMax = matchedMax
  #))
}