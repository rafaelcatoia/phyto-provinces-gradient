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

list.files(datadir)
###############################
## Loading
###############################
df_geo_abiotics = readRDS(file = paste0(savingdir,'/','df_geo_abiotics'))
dfGrump_longer_filtered = data.table::fread(paste0(datadir,'/','phyto_gradients.csv'))
list_AitDist = readRDS(paste0(savingdir,'/','list_AitDist_IS_provinces'))
list_normalized_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_abio_bio_geo_dist_normallized'))
list_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_abio_bio_geo_dist'))
grid_base = readRDS(file = paste0(savingdir,'/','grid_base'))


##
## ----------------------------------------------------------------------------- 10 neigh
##

list_cluster_membership_and_bounderies_mirroredLat = parallel::mclapply(
  list_AitDist,
  function(x){
    coloring_map(D = x,
                 nclusters = 10,
                 gbase = grid_base,
                 list_normalized_dist = list_normalized_geo_abiotics_dists)
  },mc.cores = 10)


saveRDS(list_cluster_membership_and_bounderies_mirroredLat,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_10neigh'))


##
## ----------------------------------------------------------------------------- 5 neigh
##

list_cluster_membership_and_bounderies_mirroredLat_n5 = parallel::mclapply(
  list_AitDist,
  function(x){coloring_map(D = x,nclusters = 10,
                           gbase = grid_base %>% select(lat_grid,depth_grid,n_neighs01,n_neighs02,n_neighs03,n_neighs04,n_neighs05),
                           list_normalized_dist = list_normalized_geo_abiotics_dists)},mc.cores = 10)

saveRDS(list_cluster_membership_and_bounderies_mirroredLat_n5,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_05neigh'))

##
## ----------------------------------------------------------------------------- 3 neigh
##

list_cluster_membership_and_bounderies_mirroredLat_n3 = parallel::mclapply(
  list_AitDist,
  function(x){
    coloring_map(
      D = x,nclusters = 10,
      gbase = grid_base %>% select(lat_grid,depth_grid,n_neighs01,n_neighs02,n_neighs03),
      list_normalized_dist = list_normalized_geo_abiotics_dists)},
  mc.cores = 10)

saveRDS(list_cluster_membership_and_bounderies_mirroredLat_n3,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_mirroredLat_n3'))



## Now with a different number of clusters

##
## ----------------------------------------------------------------------------- 10 neigh
##

list_cluster_membership_and_bounderies_mirroredLat = parallel::mclapply(
  list_AitDist,
  function(x){
    coloring_map(D = x,
                 nclusters = 5,
                 gbase = grid_base,
                 list_normalized_dist = list_normalized_geo_abiotics_dists)
  },mc.cores = 10)


saveRDS(list_cluster_membership_and_bounderies_mirroredLat,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_10neigh_k5'))


##
## ----------------------------------------------------------------------------- 5 neigh
##

list_cluster_membership_and_bounderies_mirroredLat_n5 = parallel::mclapply(
  list_AitDist,
  function(x){coloring_map(
    D = x,nclusters = 5,
    gbase = grid_base %>% select(lat_grid,depth_grid,n_neighs01,n_neighs02,n_neighs03,n_neighs04,n_neighs05),
    list_normalized_dist = list_normalized_geo_abiotics_dists)},mc.cores = 10)

saveRDS(list_cluster_membership_and_bounderies_mirroredLat_n5,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_05neigh_k5'))

##
## ----------------------------------------------------------------------------- 3 neigh
##

list_cluster_membership_and_bounderies_mirroredLat_n3 = parallel::mclapply(
  list_AitDist,
  function(x){
    coloring_map(
      D = x,nclusters = 5,
      gbase = grid_base %>% select(lat_grid,depth_grid,n_neighs01,n_neighs02,n_neighs03),
      list_normalized_dist = list_normalized_geo_abiotics_dists)},
  mc.cores = 10)

saveRDS(list_cluster_membership_and_bounderies_mirroredLat_n3,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_mirroredLat_n3_k5'))



