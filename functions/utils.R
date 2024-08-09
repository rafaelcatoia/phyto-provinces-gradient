## function to calculate the distance between each point in the grid and all the samples
## df_geo is the object that has scaled lat and depth
dist_grid_sample <- function(lat_depth_grid,df_geo=df_geo_abiotics,n_neigh=3){
  geo_lat_depth  = df_geo_abiotics %>% select(lat_scaled,depht_scaled) %>% as.matrix()
  vet_dist = apply(geo_lat_depth, 1, function(x){dist( rbind(x,lat_depth_grid))} )
  idx = order(vet_dist)[1:n_neigh]
  return(idx)
}


# Getting the cluster membership from df_cluster
getClustMembership <- function(x_idx,vector_cluster_mebership){
  return(vector_cluster_mebership[x_idx])
}
