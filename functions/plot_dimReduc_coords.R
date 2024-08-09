plot_dimReduc_coords <- function(coord_plots,weights = NULL,clusters=NULL,baseSize=8){
  out<-list()
  
  if(is.null(clusters)){
    if(is.null(weights)){
      out$MDS_2d = ggplot(
        data = coord_plots,
        mapping = aes(x=MDS1,y=MDS2))+
        geom_point(alpha=0.7,size=3) +
        theme_minimal(base_size = baseSize) +
        xlab(paste('MDS1')) +#,pct_explained[1],'%'))+
        ylab(paste('MDS2')) +#,pct_explained[2],'%'))+
        theme(legend.position = 'bottom')+
        ggtitle('MDS')
      
      out$nMDS_2d = ggplot(
        data = coord_plots,
        mapping = aes(x=nMDS1,y=nMDS2))+
        geom_point(alpha=0.7,size=3) +
        theme_minimal(base_size = baseSize) +
        xlab(paste('nMDS1')) +#,pct_explained[1],'%'))+
        ylab(paste('nMDS2')) +#,pct_explained[2],'%'))+
        theme(legend.position = 'bottom')+
        ggtitle('nMDS')
      
      
      out$TSNE_2d <- ggplot(
        data = coord_plots,
        mapping = aes(x=tsne1,y=tsne2))+
        geom_point(alpha=0.5,size=3) +
        theme_minimal(base_size = baseSize)+
        theme(legend.position = 'bottom')+
        ggtitle('tsne')
      
      
      out$umap_2d <- ggplot(
        data = coord_plots,
        mapping = aes(x=umap1,y=umap2))+
        geom_point(alpha=0.5,size=3) +
        theme_minimal(base_size = baseSize)+
        theme(legend.position = 'bottom')+
        ggtitle('umap')
      
      return(out)
    }else{
      coord_plots = coord_plots %>% mutate(
        weight = weights
      )
      
      out$MDS_2d = ggplot(
        data = coord_plots,
        mapping = aes(x=MDS1,y=MDS2,size=weight))+
        geom_point(alpha=0.7) +
        theme_minimal(base_size = baseSize) +
        xlab(paste('MDS1')) +#,pct_explained[1],'%'))+
        ylab(paste('MDS2')) +#,pct_explained[2],'%'))+
        theme(legend.position = 'bottom')+
        ggtitle('MDS')
      
      out$nMDS_2d = ggplot(
        data = coord_plots,
        mapping = aes(x=nMDS1,y=nMDS2,size=weight))+
        geom_point(alpha=0.7) +
        theme_minimal(base_size = baseSize) +
        xlab(paste('nMDS1')) +#,pct_explained[1],'%'))+
        ylab(paste('nMDS2')) +#,pct_explained[2],'%'))+
        theme(legend.position = 'bottom')+
        ggtitle('nMDS')
      
      
      out$TSNE_2d <- ggplot(
        data = coord_plots,
        mapping = aes(x=tsne1,y=tsne2,size=weight))+
        geom_point(alpha=0.5) +
        theme_minimal(base_size = baseSize)+
        theme(legend.position = 'bottom')+
        ggtitle('tsne')
      
      
      out$umap_2d <- ggplot(
        data = coord_plots,
        mapping = aes(x=umap1,y=umap2,size=weight))+
        geom_point(alpha=0.5) +
        theme_minimal(base_size = baseSize)+
        theme(legend.position = 'bottom')+
        ggtitle('umap')
      
      return(out)
    }
  }else{
    coord_plots = coord_plots %>% mutate(Clusters = as.factor(clusters))
    if(is.null(weights)){
      out$MDS_2d = ggplot(
        data = coord_plots,
        mapping = aes(x=MDS1,y=MDS2,color=Clusters))+
        geom_point(alpha=0.7,size=3) +
        theme_minimal(base_size = baseSize) +
        xlab(paste('MDS1')) +#,pct_explained[1],'%'))+
        ylab(paste('MDS2')) +#,pct_explained[2],'%'))+
        theme(legend.position = 'bottom')+
        ggtitle('MDS')
      
      out$nMDS_2d = ggplot(
        data = coord_plots,
        mapping = aes(x=nMDS1,y=nMDS2,color=Clusters))+
        geom_point(alpha=0.7,size=3) +
        theme_minimal(base_size = baseSize) +
        xlab(paste('nMDS1')) +#,pct_explained[1],'%'))+
        ylab(paste('nMDS2')) +#,pct_explained[2],'%'))+
        theme(legend.position = 'bottom')+
        ggtitle('nMDS')
      
      
      out$TSNE_2d <- ggplot(
        data = coord_plots,
        mapping = aes(x=tsne1,y=tsne2,color=Clusters))+
        geom_point(alpha=0.5,size=3) +
        theme_minimal(base_size = baseSize)+
        theme(legend.position = 'bottom')+
        ggtitle('tsne')
      
      
      out$umap_2d <- ggplot(
        data = coord_plots,
        mapping = aes(x=umap1,y=umap2,color=Clusters))+
        geom_point(alpha=0.5,size=3) +
        theme_minimal(base_size = baseSize)+
        theme(legend.position = 'bottom')+
        ggtitle('umap')
      
      return(out)
    }else{
      coord_plots = coord_plots %>% mutate(
        weight = weights
      )
      
      out$MDS_2d = ggplot(
        data = coord_plots,
        mapping = aes(x=MDS1,y=MDS2,size=weight,color=Clusters))+
        geom_point(alpha=0.7) +
        theme_minimal(base_size = baseSize) +
        xlab(paste('MDS1')) +#,pct_explained[1],'%'))+
        ylab(paste('MDS2')) +#,pct_explained[2],'%'))+
        theme(legend.position = 'bottom')+
        ggtitle('MDS')
      
      out$nMDS_2d = ggplot(
        data = coord_plots,
        mapping = aes(x=nMDS1,y=nMDS2,size=weight,color=Clusters))+
        geom_point(alpha=0.7) +
        theme_minimal(base_size = baseSize) +
        xlab(paste('nMDS1')) +#,pct_explained[1],'%'))+
        ylab(paste('nMDS2')) +#,pct_explained[2],'%'))+
        theme(legend.position = 'bottom')+
        ggtitle('nMDS')
      
      
      out$TSNE_2d <- ggplot(
        data = coord_plots,
        mapping = aes(x=tsne1,y=tsne2,size=weight,color=Clusters))+
        geom_point(alpha=0.5) +
        theme_minimal(base_size = baseSize)+
        theme(legend.position = 'bottom')+
        ggtitle('tsne')
      
      
      out$umap_2d <- ggplot(
        data = coord_plots,
        mapping = aes(x=umap1,y=umap2,size=weight,color=Clusters))+
        geom_point(alpha=0.5) +
        theme_minimal(base_size = baseSize)+
        theme(legend.position = 'bottom')+
        ggtitle('umap')
      
      return(out)
    }
  }
}

