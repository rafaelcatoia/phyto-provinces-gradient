####################################################### ---
library(parallel) ; library(dplyr) ; library(tidyr)

root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")

files_vec <- list.files(funsdir)

for( i in 1:length(files_vec)){
  source(root$find_file(paste0(funsdir,'/',files_vec[i])))
}

df_geo_abiotics <- readRDS(paste0(savingdir,'/','df_geo_abiotics'))
df_geo_abiotics = df_geo_abiotics %>%
  mutate(SampleID=factor(SampleID)) %>% 
  arrange(SampleID)

## Store mean and sd of the df_geo_abiotics
lat_scale_obj   <- c(mean(df_geo_abiotics$Latitude),sd(df_geo_abiotics$Latitude))
depth_scale_obj <- c(mean(df_geo_abiotics$Depth),sd(df_geo_abiotics$Depth))

df_geo_abiotics = df_geo_abiotics %>% 
  mutate(lat_scaled = (Latitude - mean(Latitude) ) /sd(Latitude),
         depht_scaled = (Depth -mean(Depth) ) / sd(Depth))

## Now setting the values that the grid will run on
min_lat = round(min(df_geo_abiotics$lat_scaled),1)-0.1
max_lat = round(max(df_geo_abiotics$lat_scaled),1)+0.1
min_depth = round(min(df_geo_abiotics$depht_scaled),1)-0.1
max_depth = round(max(df_geo_abiotics$depht_scaled),1)+0.1

lat_grid = seq(min_lat,max_lat,0.05)
depth_grid = seq(min_depth,max_depth,0.05)

## expanding the grid 
grid_base = expand_grid(lat_grid,depth_grid)

## Here we set the number of nearest neigh 
nneigh = 10

#matrix that will retain the index of the neighbors
nneigh_aux <- matrix(NA,nrow = nrow(grid_base),ncol=nneigh)
maxIt = nrow(grid_base)

for( i in 1:maxIt){
  nneigh_aux[i,] <- dist_grid_sample(grid_base[i,],n_neigh = nneigh)
  if(i%%100==0){cat('iteration ----------------------- ',i,'of ',maxIt,'\n')}
}

## only giving names to the columns 
aux <- 1:ncol(nneigh_aux)
colnames(nneigh_aux) = ifelse(aux<10,
                              paste0('n_neighs0',aux),
                              paste0('n_neighs',aux))


grid_base = bind_cols(grid_base %>% select(lat_grid,depth_grid),nneigh_aux) 

## now we can basically start from here now.
saveRDS(grid_base,paste0(savingdir,'/','grid_base'))
#grid_base <- readRDS(paste0(savingdir,'/','grid_base'))

## Everything underneath was only used to create the functions required to cloloring/finding the limits.

# ###### Get sample cluster membership based on index
# nneigh_colnames <- grid_base %>% select(starts_with('n_neig')) %>% colnames()
# aux = 1:length(nneigh_colnames)
# custDominates =  ifelse(aux<10,
#                               paste0('ClustDom0',aux),
#                               paste0('ClustDom',aux))
# matCluster = matrix(NA,nrow=nrow(grid_base),ncol=length(nneigh_colnames))
# 
# 
# for(i in 1:length(nneigh_colnames)){
#   clustDom = grid_base %>% select(any_of(nneigh_colnames[i])) %>% pull()
#   matCluster[,i] <- getClustMembership(x_idx = clustDom,vector_cluster_mebership = df_cluster$Cluster)
# }
# 
# 
# ## here we have a possible problem. Some points have more than 
# grid_base = grid_base %>% mutate(ClustRegion = most_frequent_value(matCluster))
# 
# grid_base_matrix = grid_base %>% 
#   pivot_wider(id_cols = lat_grid, names_from = depth_grid, values_from = ClustRegion) %>% 
#   as.matrix()
# 
# grid_base_matrix_limits = check_adjacent_values(grid_base_matrix[,-1]) %>% data.frame() %>% 
#   pivot_longer(everything())
# 
# grid_base = grid_base %>% 
#   mutate(Lat=lat_grid*lat_scale_obj[2] + lat_scale_obj[1]   ,
#          Depth = depth_grid*depth_scale_obj[2] + depth_scale_obj[1],
#          limits= rowMeans(df_Limits %>% filter(clust_control=="ward_0") %>% select(-clust_control))
#          )
# 
# colourCount = length(df_cluster$Cluster %>% unique())
# getPalette = colorRampPalette(brewer.pal(9, "Set1"))
# 
# pal1 = RColorBrewer::brewer.pal(colourCount,)
# 
# mypalette<-RColorBrewer::brewer.pal(12,"Paired")
# 
# ggplot()+
#   #geom_tile(aes(x=Lat,y=Depth,fill = factor(ClustRegion)),
#   #          data = grid_base %>% filter(Depth>=0),alpha=0.5)+
#   geom_tile(data = grid_base %>% filter(Depth>=0) %>% filter(limits>0),
#              aes(x=Lat,y=Depth,fill=limits,alpha=limits),size=0.1)+
#   geom_label(aes(x=Latitude,y=Depth,label=Cluster,color=Cluster),
#              data = df_cluster %>% mutate(Cluster=as.factor(Cluster)) %>% 
#                left_join(df_geo_abiotics %>% transmute(SampleID,Latitude,Depth)))+
#   theme_minimal()+
#   scale_y_reverse()+
#   theme(legend.position = 'bottom')
# 
# grid_base
# 
# ggplot()+
#   geom_tile(aes(x=lat_grid,y=depth_grid,fill = factor(ClustRegion)),
#             data = grid_base %>% filter(Depth>=0),alpha=0.5)+
#   geom_point(data = grid_base %>% filter(Depth>=0) %>% filter(!is.na(limits)),
#              aes(x=lat_grid,y=depth_grid,fill=factor(limits)),size=0.1,color='white')+
#   geom_label(aes(x=Latitude,y=Depth,label=Cluster,color=Cluster),
#              data = df_cluster %>% mutate(Cluster=as.factor(Cluster)) %>% 
#                left_join(df_geo_abiotics %>% 
#                            transmute(SampleID,
#                                      Latitude=(Latitude-mean(Latitude))/sd(Latitude),
#                                      Depth=(Depth-mean(Depth))/sd(Depth))))+
#   theme_minimal()+
#   scale_y_reverse()+
#   theme(legend.position = '')
# 
# 
# ########## -------------------------------------------------------------
# p1 = ggplot()+
#   geom_tile(aes(x=Lat,y=Depth,fill = factor(ClustRegion)),
#             data = grid_base %>% filter(Depth>=0),alpha=0.5)+
#   theme_minimal()+
#   scale_y_reverse()+
#   theme(legend.position = 'bottom')
# 
# p1
# 
# p2 = p1 +  geom_point(data = grid_base %>% filter(Depth>=0) %>% filter(limits>0.25),
#                       aes(x=Lat,y=Depth,color=limits,alpha=limits),size=0.25)+
#   scale_color_gradient(low ='white',high = 'white')
# 
# p2
# 
# p3 = p2 +ggnewscale::new_scale(new_aes = 'Cluster')  +
#   
#   geom_label(
#     aes(x=Latitude,y=Depth,label=Cluster,color=Cluster),
#     data = df_cluster %>% mutate(Cluster=as.factor(Cluster)) %>% 
#       left_join(df_geo_abiotics %>% 
#                   transmute(SampleID,
#                             Latitude ,#=#(Latitude-mean(Latitude))/sd(Latitude),
#                             Depth#=(Depth-mean(Depth))/sd(Depth))))
#                   )))+
#   scale_colour_discrete() +
#   labs(colour = "Black-Grey",y = "Cluster")
