################################## 
library(dplyr) ; library(tidyr)

root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")

files_vec <- list.files(funsdir)

for( i in 1:length(files_vec)){
  source(root$find_file(paste0(funsdir,'/',files_vec[i])))
}

###############################
## Loading
###############################
df_geo_abiotics = readRDS(file = paste0(savingdir,'/','df_geo_abiotics'))
#dfGrump_longer_filtered = data.table::fread(paste0(datadir,'/','grump_phyto.csv'))
list_AitDist = readRDS(paste0(savingdir,'/','list_AitDist_IS_provinces'))
list_normalized_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_abio_bio_geo_dist_normallized'))
list_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_abio_bio_geo_dist'))
grid_base = readRDS(file = paste0(savingdir,'/','grid_base'))

## Creating cluster label ####### 
nclusters_match = 3
D = list_normalized_geo_abiotics_dists$bioticDist
vet_alpha = c(0,0.05,0.1,0.2)
D1 = as.dist(D)
D2 = as.dist( (1-as.numeric(vet_alpha[1]))*D + as.numeric(vet_alpha[1])*list_normalized_geo_abiotics_dists$geoDist)
D3 = as.dist( (1-as.numeric(vet_alpha[2]))*D + as.numeric(vet_alpha[2])*list_normalized_geo_abiotics_dists$geoDist)
D4 = as.dist( (1-as.numeric(vet_alpha[3]))*D + as.numeric(vet_alpha[3])*list_normalized_geo_abiotics_dists$geoDist)

####################################   Ward    ####################################
#################################### Creating hclust objs
hclust_0 = hclust(d=D1,method='ward.D')
hclust_0.05 = hclust(d=D2,method='ward.D')
hclust_0.10 = hclust(d=D3,method='ward.D')
hclust_0.20 = hclust(d=D4,method='ward.D')

#################################### With geodist 
mat_cluster_membership_label <- cbind.data.frame(
  v1_hclust=cutree(hclust_0,k = nclusters),
  v2_hclust=cutree(hclust_0.05,k = nclusters),
  v3_hclust=cutree(hclust_0.10,k = nclusters),
  v4_hclust=cutree(hclust_0.20,k = nclusters)
)
######## ---------------------------------------------------------------------
## Getting only part of the chain
list_AitDist <- list_AitDist[1:100]


## Lets do one for sanity checking
# testing <- coloring_map_matching(
#   D=list_AitDist[[1]],
#   nclusters = 3,
#   trueClusterMembership = mat_cluster_membership_label,
#   gbase = grid_base,
#   list_normalized_dist = list_normalized_geo_abiotics_dists)


## ----------------------------------------------------------------------------- 10 neigh
##

list_cluster_membership_and_bounderies_mirroredLat = parallel::mclapply(
  list_AitDist,
  function(x){
    coloring_map_matching(
      D = x,
      nclusters = nclusters_match,
      trueClusterMembership = mat_cluster_membership_label,
      gbase = grid_base,
      list_normalized_dist = list_normalized_geo_abiotics_dists)
  },mc.cores = 10)


saveRDS(list_cluster_membership_and_bounderies_mirroredLat,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_10neigh_matched'))


##
## ----------------------------------------------------------------------------- 5 neigh
##

list_cluster_membership_and_bounderies_mirroredLat_n5 = parallel::mclapply(
  list_AitDist,
  function(x){
    coloring_map_matching(
      D = x,nclusters = nclusters_match,
      trueClusterMembership = mat_cluster_membership_label,
      gbase = grid_base %>% select(lat_grid,depth_grid,n_neighs01,n_neighs02,n_neighs03,n_neighs04,n_neighs05),
      list_normalized_dist = list_normalized_geo_abiotics_dists)
  },
  mc.cores = 10)

saveRDS(list_cluster_membership_and_bounderies_mirroredLat_n5,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_05neigh_matched'))

##
## ----------------------------------------------------------------------------- 3 neigh
##

list_cluster_membership_and_bounderies_mirroredLat_n3 = parallel::mclapply(
  list_AitDist,
  function(x){
    coloring_map_matching(
      D = x,nclusters = nclusters_match,
      trueClusterMembership = mat_cluster_membership_label,
      gbase = grid_base %>% select(lat_grid,depth_grid,n_neighs01,n_neighs02,n_neighs03),
      list_normalized_dist = list_normalized_geo_abiotics_dists)},
  mc.cores = 10)

saveRDS(list_cluster_membership_and_bounderies_mirroredLat_n3,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_mirroredLat_n3_matched'))



