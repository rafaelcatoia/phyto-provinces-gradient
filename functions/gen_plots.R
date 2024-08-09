######## Genete plots ---------------------------------
# list_cluster_membership_and_bounderies = readRDS(file = paste0(savingdir,'/','list_cluster_membership_and_bounderies'))
# list_cluster_membership_and_bounderies_mirroredLat = readRDS(file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_mirroredLat'))
# grid_base = readRDS(file = paste0(savingdir,'/','grid_base'))
# clust_member_limits = unlist_coloring_obj(list_cluster_membership_and_bounderies_mirroredLat)
# df_GeoAbio = readRDS(file = paste0(savingdir,'/','df_geo_abiotics'))

#gridBase = grid_base %>% rename('lat_grid'=1,'depth_grid'=2)
#clustMemLim = clust_member_limits
#df_GeoAbio =df_geo_abiotics

gen_plots <- function(gridBase,clustMemLim,df_GeoAbio,
                      LongHurst = F,oceanNames=F,verticalSlices=F,
                      plotClusterMembership = NULL){
  
  vet_longhurst <- c(-67.5, -57.2, - 44, -38, -7.5, 4 ,10, 34.5, 44.5)
  vet_oceanSlices <- c(-60,-40,-5,5,35)
  if(LongHurst){
  vet_vlines <- vet_longhurst
  }else{
    vet_vlines <- vet_oceanSlices
  }
  
  ## first lets bring LAT and Depht to it's original scale
  mean_Lat = df_GeoAbio$Latitude %>% mean()
  sd_Lat = df_GeoAbio$Latitude %>% sd()
  
  mean_Depth = df_GeoAbio$Depth %>% mean()
  sd_Depth = df_GeoAbio$Depth %>% sd()

  gridBase = gridBase %>% 
    mutate(Latitude= lat_grid*sd_Lat + mean_Lat,
           Depth = depth_grid*sd_Depth + mean_Depth)
  
  #those values are going to use as limits on the plots
  min_lat = min(gridBase$Latitude)
  max_lat = max(gridBase$Latitude)
  min_Depth = min(gridBase$Depth)
  max_Depth = max(gridBase$Depth)
  #################### First the summary plots ----------------------------------------------
  meanLimits = lapply(clustMemLim$regionLimits,rowMeans)
  clustRegion = parallel::mclapply(clustMemLim$regionCluster,most_frequent_random,mc.cores = 10)
  number_replicates = ncol(clustMemLim$regionCluster[[1]])
  labelSize=12
  
  #### Plotting limits 
  list_summary_Limits = lapply(seq_along(meanLimits), function(i){
    
    plt = gridBase %>% select(Latitude,Depth) %>% 
      mutate(limits = meanLimits[[i]]) %>% filter(limits>0) %>% 
      filter(Depth > 0 ) %>% 
      ggplot()+
      geom_tile(aes(x=Latitude,y=Depth,fill=limits))+
      scale_y_reverse()+
      theme_minimal()+
      ggtitle(names(meanLimits)[i])+
      xlim(min_lat,max_lat)+
      scale_fill_gradient(low = "white", high = "black")+
      theme(legend.position ="bottom")
      #theme(legend="bottom")
    
    if(verticalSlices){
      plt = plt + geom_vline(xintercept = vet_vlines,col='deepskyblue3',alpha=0.7)
    }
    
    if(oceanNames){
      plt = plt +
      annotate("text",size = labelSize*0.3, x = -68, y = -25, label = "Southern Ocean",hjust=0.5)+
      annotate("text",size = labelSize*0.3, x = -50, y = -25, label = "Subantarctic",hjust=0.5)+
      annotate("text",size = labelSize*0.3, x = -22, y = -25, label = "South Pacific Gyre",hjust=0.5)+
      annotate("text",size = labelSize*0.3, x = 0,   y = -25, label = "Equatorial",hjust=0.5)+
      annotate("text",size = labelSize*0.3, x = 20,  y = -25, label = "North Pacific Gyre",hjust=0.5)+
      annotate("text",size = labelSize*0.3, x = 50,  y = -25, label = "Subarctic",hjust=0.5)
    }
    
    if(!is.null(plotClusterMembership)){
      df_plot_clust = data.frame(plotClusterMembership)[,c(1,2,i+2)]
      colnames(df_plot_clust)<- c("Latitude","Depth","Cluster")
      df_plot_clust = df_plot_clust %>% mutate(Cluster=as.factor(Cluster))
      plt <- plt + 
        geom_text(data = df_plot_clust,aes(x=Latitude,y=Depth,label=Cluster))
    }
    
    return(plt)
  })
  
  limits_faceted <- ggpubr::ggarrange(plotlist = list_summary_Limits,ncol = 4,nrow=2,
                                      common.legend = T, legend="bottom")
  
  
  list_summary_clustRegion = lapply(seq_along(clustRegion), function(i){
    
    plt = gridBase %>% select(Latitude,Depth) %>% 
      mutate(clust_Region = factor(clustRegion[[i]]$value),
             pct = clustRegion[[i]]$frequency/number_replicates) %>% 
      filter(Depth >0 ) %>% 
      ggplot()+
      geom_tile(aes(x=Latitude,y=Depth,fill=clust_Region,alpha=pct))+
      scale_y_reverse()+
      theme_minimal()+
      ggtitle(names(clustRegion)[i]) %>% 
      theme(legend.position = 'bottom')
    
    if(verticalSlices){
        plt = plt + geom_vline(xintercept = vet_vlines,col='deepskyblue3')
    }
    
    if(oceanNames){
      plt = plt +
        annotate("text",size = labelSize*0.3, x = -68, y = -25, label = "Southern Ocean",hjust=0.5)+
        annotate("text",size = labelSize*0.3, x = -50, y = -25, label = "Subantarctic",hjust=0.5)+
        annotate("text",size = labelSize*0.3, x = -22, y = -25, label = "South Pacific Gyre",hjust=0.5)+
        annotate("text",size = labelSize*0.3, x = 0,   y = -25, label = "Equatorial",hjust=0.5)+
        annotate("text",size = labelSize*0.3, x = 20,  y = -25, label = "North Pacific Gyre",hjust=0.5)+
        annotate("text",size = labelSize*0.3, x = 50,  y = -25, label = "Subarctic",hjust=0.5)}
    
    if(!is.null(plotClusterMembership)){
      df_plot_clust = data.frame(plotClusterMembership)[,c(1,2,i+2)]
      colnames(df_plot_clust)<- c("Latitude","Depth","Cluster")
      df_plot_clust = df_plot_clust %>% mutate(Cluster=as.factor(Cluster))
      plt <- plt + 
        geom_text(data = df_plot_clust,aes(x=Latitude,y=Depth,label=Cluster))
    }
    
    return(plt)
    
    
    })
  
  ## Had to do this in order to have the same 
  legend_unique = ggpubr::get_legend(list_summary_clustRegion[[4]])
  
  clustRegion_facetd <- ggpubr::ggarrange(
    plotlist = list_summary_clustRegion,common.legend = T,legend.grob = legend_unique,
    ncol = 4,nrow=2, legend="bottom")
  
  return(list(
    list_summary_clustRegion = list_summary_clustRegion,
    list_limits = list_summary_Limits,
    limits_faceted=limits_faceted,
    clustRegion_facetd=clustRegion_facetd
  )
  )  
  
}

most_frequent_random <- function(matrix) {
  get_most_frequent <- function(row) {
    tab <- table(row)
    max_count <- max(tab)
    most_frequent_values <- names(tab[tab == max_count])
    selected_value <- sample(most_frequent_values, 1)
    return(c(as.numeric(selected_value), max_count))
  }
  
  result <- t(apply(matrix, 1, get_most_frequent))
  df <- data.frame(value = result[, 1], frequency = result[, 2])
  return(df)
}

## Example usage
#matrix <- matrix(c(1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5), nrow = 4, byrow = TRUE)
#print(matrix)
#print(most_frequent(matrix))