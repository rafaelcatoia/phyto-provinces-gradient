##### Matching Function 

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

relabel_matching <- function(x,xMatched){
  out_x <- numeric(length = length(x))
  for(i in 1:max(xMatched)){
    out_x[which(x==xMatched[i])] <- i
  }
  return(out_x)
}


# batata <- function(matX){
#   dimMat = dim(as.matrix(matX))
#   #if no dimension is 1 
#   if(!(sum(dimMat==1)>0)){
#     if(dimMat[1]>dimMat[2]){
#       return(mean(apply(matX,1,min)))
#     }else{
#       return(mean(apply(matX,2,min)))
#     }
#   }else{
#     return(mean(as.vector(matX)))
#   }
# }

## Checking ----
# 
# obj_matched = matching_clusters(
#   DistMatrix = list_geo_abiotics_dists$geo_Dist_mirrored,
#   trueLabel = trueClusterMembership,
#   proposedLabel = current_label
# )
# mean(current_label==trueClusterMembership)
# mean(obj_matched$matchedMean==trueClusterMembership)
# 
# p1 = df_geo_abiotics %>% mutate(
#   trueLabel = trueClusterMembership,
#   replicateLabel = as.factor(current_label)
# ) %>% ggplot(
#   aes(x=Latitude,y=Depth)
# )+geom_label(aes(label=trueLabel,color=replicateLabel))+
#   scale_y_reverse()
# 
# p2 = df_geo_abiotics %>% mutate(
#   trueLabel = trueClusterMembership,
#   replicateLabel = as.factor(obj_matched$matchedMean)
# ) %>% ggplot(
#   aes(x=Latitude,y=Depth)
# )+geom_label(aes(label=trueLabel,color=replicateLabel))+
#   scale_y_reverse()
# 
# p3 = df_geo_abiotics %>% mutate(
#   trueLabel = trueClusterMembership,
#   replicateLabel = as.factor(obj_matched$matchedMax)
# ) %>% ggplot(
#   aes(x=Latitude,y=Depth)
# )+geom_label(aes(label=trueLabel,color=replicateLabel))+
#   scale_y_reverse()
# 
# ggpubr::ggarrange(plotlist = list(p1,p2,p3),
#                   nrow = 1,
#                   common.legend = T)
# 