############# ----------------------------------------------------------
## script to create a function to evaluate the clusterings -- deciding K
############# ----------------------------------------------------------
 
# df_geo_abiotics = readRDS(file = paste0(savingdir,'/','df_geo_abiotics'))
# dfGrump_longer_filtered = data.table::fread(paste0(datadir,'/','grump_phyto.csv'))
# list_AitDist = readRDS(file = paste0(savingdir,'/','list_AitDist_IS_OOS'))
# list_normalized_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_abio_bio_geo_dist_normallized'))
# list_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_abio_bio_geo_dist'))
# D = list_AitDist[[1]]

eval_clustering_OOS <- function(
  D,
  list_normalized_dist=list_normalized_geo_abiotics_dists, # those dists can be used to create a convex mixture of distance matrices
  list_dist = list_geo_abiotics_dists, # those will be use to perform cluser evaluation
  vet_alpha=c(0.05,0.10,0.20),
  Kmax=30){
  
  #D$INS=as.matrix(D$INS)
  D_matrix = as.matrix(D$INS)
  
  nrows = dim(D_matrix)[1]
  ####################################                ####################################
  #################################### Creating Phase ####################################
  ####################################                ####################################
  
  ### Settind dist matrices 
  

  ### Settind dist matrices 
  D1 = as.dist(D_matrix)
  D2 = as.dist( (1-as.numeric(vet_alpha[1]))*D_matrix + as.numeric(vet_alpha[1])*list_normalized_dist$geoDist)
  D3 = as.dist( (1-as.numeric(vet_alpha[2]))*D_matrix + as.numeric(vet_alpha[2])*list_normalized_dist$geoDist)
  D4 = as.dist( (1-as.numeric(vet_alpha[3]))*D_matrix + as.numeric(vet_alpha[3])*list_normalized_dist$geoDist)
  
  ####################################   Ward    ####################################
  #################################### Creating hclust objs
  hclust_0 = hclust(d=D1,method='ward.D')
  hclust_0.05 = hclust(d=D2,method='ward.D')
  hclust_0.10 = hclust(d=D3,method='ward.D')
  hclust_0.20 = hclust(d=D4,method='ward.D')
  
  #################################### With geodist 
  mat_hclust0 <- matrix(NA,nrow=nrows,ncol=Kmax)
  mat_hclust0[,1]<-1
  
  mat_hclust_0.05<- matrix(NA,nrow=nrows,ncol=Kmax)
  mat_hclust_0.05[,1]<-1
  
  mat_hclust_0.1<- matrix(NA,nrow=nrows,ncol=Kmax)
  mat_hclust_0.1[,1]<-1
  
  mat_hclust_0.2<- matrix(NA,nrow=nrows,ncol=Kmax)
  mat_hclust_0.2[,1]<-1
  
  for(ii in 2:Kmax){
    mat_hclust0[,ii]=cutree(hclust_0,k = ii)
    mat_hclust_0.05[,ii]=cutree(hclust_0.05,k = ii)
    mat_hclust_0.1[,ii]=cutree(hclust_0.10,k = ii)
    mat_hclust_0.2[,ii]=cutree(hclust_0.20,k = ii)
  }
  
  ## We are not doing tha for pam ---
  # ####################################   pam    ####################################
  # mat_pam0 <- matrix(NA,nrow=nrows,ncol=Kmax)
  # mat_pam0[,1]<-1
  # 
  # mat_pam_0.05<- matrix(NA,nrow=nrows,ncol=Kmax)
  # mat_pam_0.05[,1]<-1
  # 
  # mat_pam_0.1<- matrix(NA,nrow=nrows,ncol=Kmax)
  # mat_pam_0.1[,1]<-1
  # 
  # mat_pam_0.2<- matrix(NA,nrow=nrows,ncol=Kmax)
  # mat_pam_0.2[,1]<-1

  
  # for(ii in 2:Kmax){
  #   mat_pam0[,ii]     = cluster::pam(x=D1,k = ii,cluster.only = T)
  #   mat_pam_0.05[,ii] = cluster::pam(x=D2,k = ii,cluster.only = T)
  #   mat_pam_0.1[,ii]  = cluster::pam(x=D3,k = ii,cluster.only = T)
  #   mat_pam_0.2[,ii]  = cluster::pam(x=D4,k = ii,cluster.only = T)
  # }
  
  
  ####################################                  ####################################
  #################################### Evaluating Phase ####################################
  ####################################                  ####################################
  Kmax = max(mat_hclust0)
  
  methodId = 'ward'
  df_hclust = bind_rows(
    bind_rows(
      evaluate_cluster(distMatrixEval = list_dist$geoDist %>% as.matrix,df_clusters = mat_hclust0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
      evaluate_cluster(distMatrixEval = list_dist$geoDist %>% as.matrix,df_clusters = mat_hclust_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
      evaluate_cluster(distMatrixEval = list_dist$geoDist %>% as.matrix,df_clusters = mat_hclust_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
      evaluate_cluster(distMatrixEval = list_dist$geoDist %>% as.matrix,df_clusters = mat_hclust_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
    ) %>% mutate(DistMetric='Geo') ,
    bind_rows(
      evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_hclust0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
      evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_hclust_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
      evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_hclust_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
      evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_hclust_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
    ) %>% mutate(DistMetric='Abio'),
    bind_rows(
      evaluate_cluster(distMatrixEval = D$OOS ,df_clusters = mat_hclust0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
      evaluate_cluster(distMatrixEval = D$OOS ,df_clusters = mat_hclust_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
      evaluate_cluster(distMatrixEval = D$OOS ,df_clusters = mat_hclust_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
      evaluate_cluster(distMatrixEval = D$OOS ,df_clusters = mat_hclust_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
    ) %>% mutate(DistMetric='Bio')
  )
  
  # methodId ='pam'
  # df_pam = bind_rows(
  #   bind_rows(
  #     evaluate_cluster(distMatrixEval = list_dist$geo_Dist %>% as.matrix,df_clusters = mat_pam0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
  #     evaluate_cluster(distMatrixEval = list_dist$geo_Dist %>% as.matrix,df_clusters = mat_pam_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
  #     evaluate_cluster(distMatrixEval = list_dist$geo_Dist %>% as.matrix,df_clusters = mat_pam_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
  #     evaluate_cluster(distMatrixEval = list_dist$geo_Dist %>% as.matrix,df_clusters = mat_pam_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
  #   ) %>% mutate(DistMetric='Geo'),
  #   bind_rows(
  #     evaluate_cluster(distMatrixEval = list_dist$geo_Dist_mirrored %>% as.matrix,df_clusters = mat_pam0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
  #     evaluate_cluster(distMatrixEval = list_dist$geo_Dist_mirrored %>% as.matrix,df_clusters = mat_pam_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
  #     evaluate_cluster(distMatrixEval = list_dist$geo_Dist_mirrored %>% as.matrix,df_clusters = mat_pam_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
  #     evaluate_cluster(distMatrixEval = list_dist$geo_Dist_mirrored %>% as.matrix,df_clusters = mat_pam_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
  #   ) %>% mutate(DistMetric='Geo_Mirrored'),
  #   bind_rows(
  #     evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_pam0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
  #     evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_pam_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
  #     evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_pam_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
  #     evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_pam_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
  #   ) %>% mutate(DistMetric='Abiotics'),
  #   bind_rows(
  #     evaluate_cluster(distMatrixEval = D$OOS ,df_clusters = mat_pam0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
  #     evaluate_cluster(distMatrixEval = D$OOS ,df_clusters = mat_pam_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
  #     evaluate_cluster(distMatrixEval = D$OOS ,df_clusters = mat_pam_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
  #     evaluate_cluster(distMatrixEval = D$OOS ,df_clusters = mat_pam_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
  #   ) %>% mutate(DistMetric='Aitchison')
  # )
  
  #df_out = bind_rows(df_hclust,df_pam)
  df_out = df_hclust
  return(df_out) 
}


#eval_clustering_OOS(D = D,
#                    list_normalized_dist = list_normalized_geo_abiotics_dists,
#                    list_dist = list_geo_abiotics_dists,Kmax = 30,
#                    vet_alpha = c(0.05,0.1,0.2))
