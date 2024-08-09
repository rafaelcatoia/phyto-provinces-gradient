################ ----------------------------------------------------------------------
## Function to work with the output from lapply(list_AitDist_clean[1:100],coloring_map)
################ ----------------------------------------------------------------------

# list_cluster_membership_and_bounderies_mirroredLat = readRDS(file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_mirroredLat'))
# 
# list_boundries = list_cluster_membership_and_bounderies_mirroredLat

unlist_coloring_obj <- function(list_boundries){
  list_list_2_list_df = lapply(list_boundries, stack_dataframes)
  list_list_2_list_df = lapply(list_list_2_list_df,data.frame)
  
  ## now lets give those replicates some names
  names(list_list_2_list_df) <- paste0('repli_',1:length(list_list_2_list_df))
  
  ## getting the lenght of replicates
  leng_list = length(list_list_2_list_df)
  nrows = nrow(list_list_2_list_df[[1]])
  
  ## creating the outputs
  df_clustRegion <- matrix(NA,nrow =nrows ,ncol=leng_list)
  df_Limits <- matrix(NA,nrow =nrows ,ncol=leng_list)
  
  for(ii in 1:leng_list){
    df_clustRegion[,ii] = list_list_2_list_df[[ii]]$clustRegion
    df_Limits[,ii] = list_list_2_list_df[[ii]]$limits
  }
  
  df_clustRegion = cbind.data.frame(
    clust_control = list_list_2_list_df$repli_1$index_name,
    df_clustRegion
  )
  df_Limits = cbind.data.frame(
    clust_control = list_list_2_list_df$repli_1$index_name,
    df_Limits
  )
  
  names_list <- df_clustRegion$clust_control %>% unique
  names_list =  names_list[order(names_list)]
  
  list_Limits = df_Limits %>% group_split(clust_control)
  list_Limits = lapply(list_Limits,function(x){x %>% select(-clust_control)})
  names(list_Limits)<- names_list
  
  list_regions = df_clustRegion %>% group_split(clust_control)
  list_regions = lapply(list_regions,function(x){x %>% select(-clust_control)})
  names(list_regions)<- names_list
  
  ## the order of the output is pretty inportant and should be binded to a grid_base.
  return(
    list(
      regionCluster = list_regions,
      regionLimits = list_Limits
    )
  )
  
}






