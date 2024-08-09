library(dplyr) ; library(tidyr)

root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")

library(parallel) ; library(dplyr)

files_vec <- list.files(funsdir)
for( i in 1:length(files_vec)){
  source(root$find_file(paste0(funsdir,'/',files_vec[i])))
}

###############################
## Loading
###############################
df_geo_abiotics = readRDS(file = paste0(savingdir,'/','df_geo_abiotics'))
dfGrump_longer_filtered = data.table::fread(paste0(datadir,'/','phyto_gradients.csv'))
list_AitDist = readRDS(file = paste0(savingdir,'/','list_AitDist_IS_OOS'))
list_normalized_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_abio_bio_geo_dist_normallized'))
list_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_abio_bio_geo_dist'))

########################### Deciding K ###########################
## Computing the metrics ##
df_evaluation <- data.table::rbindlist(
  parallel::mclapply(
    list_AitDist,function(x){
      eval_clustering_OOS(
        D = x,
        list_normalized_dist = list_normalized_geo_abiotics_dists,
        list_dist = list_geo_abiotics_dists,
        Kmax = 50)},mc.cores = 10))

saveRDS(df_evaluation,paste0(savingdir,'/','df_evaluation_choosing_k'))
#df_evaluation = readRDS(paste0(savingdir,'/','df_evaluation'))

######### looking at it now ::

## First lets arrange all the metrics and summarise
df_summary = bind_rows(
  df_evaluation %>% 
    group_by(ncluster,method,alpha,DistMetric) %>% 
    summarise_all(.funs = mean) %>% mutate(summary_metric='Mean') %>% ungroup() %>% 
    pivot_longer(cols = -c(ncluster,method,alpha,DistMetric,summary_metric),names_to = 'Metric'),
  
  df_evaluation %>% 
    group_by(ncluster,method,alpha,DistMetric) %>% 
    summarise_all(.funs = sd) %>% mutate(summary_metric='SD') %>% ungroup() %>% 
    pivot_longer(cols = -c(ncluster,method,alpha,DistMetric,summary_metric),names_to = 'Metric')
)

library(ggplot2)

## now lets plot -----

## First Sratio
df_summary %>% filter(Metric=='avg_within_between_dist2medoid') %>% select(-Metric) %>% 
  pivot_wider(id_cols = ncluster:DistMetric,names_from = summary_metric ) %>%
  mutate(alpha=factor(alpha)) %>% 
  ggplot(aes(x=ncluster,y=Mean,color=alpha,fill=alpha,linetype=alpha))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_grid(DistMetric~alpha,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')+
  ggtitle('Sratio')


## Need Ideas
df_summary %>% filter(Metric=='avg_within_between_dist2medoid') %>%
  select(-Metric) %>% 
  pivot_wider(id_cols = ncluster:DistMetric,names_from = summary_metric ) %>%
  mutate(alpha=factor(alpha)) %>% 
  filter(alpha=='0',DistMetric=='Bio') %>% 
  mutate(CumSum = cumsum(Mean),
         Test = Mean/CumSum) %>% 
  ggplot(aes(x=ncluster,y=Test,color=alpha,fill=alpha,linetype=alpha))+
  geom_line()+
  #geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_wrap(~DistMetric,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')+
  ggtitle('Sratio')



df_summary %>% filter(Metric=='avg_dist_within_dist2medoid') %>% select(-Metric) %>% 
  pivot_wider(id_cols = ncluster:DistMetric,names_from = summary_metric ) %>%
  mutate(alpha=factor(alpha)) %>% 
  ggplot(aes(x=ncluster,y=Mean,color=alpha,fill=alpha,linetype=alpha))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_grid(DistMetric~alpha,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')+
  ggtitle('Snumerator')


df_summary %>% filter(Metric=='avg_dist_between_medoids') %>% select(-Metric) %>% 
  pivot_wider(id_cols = ncluster:DistMetric,names_from = summary_metric ) %>%
  mutate(alpha=factor(alpha)) %>% 
  ggplot(aes(x=ncluster,y=Mean,color=alpha,fill=alpha,linetype=alpha))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_grid(DistMetric~alpha,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')+
  ggtitle('Sdenon')

df_summary %>% filter(Metric=='max_max_dist_within') %>% select(-Metric) %>% 
  pivot_wider(id_cols = ncluster:DistMetric,names_from = summary_metric ) %>%
  mutate(alpha=factor(alpha)) %>% 
  ggplot(aes(x=ncluster,y=Mean,color=method,fill=method,linetype=alpha))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_grid(DistMetric~alpha,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')+
  ggtitle('max_max_dist_within')


########################### Deciding alpha ###########################
### I'll keep using K=10 for now.
### Let's focus in tuning alpha -------------------------------------------------------

df_evaluation <- data.table::rbindlist(
  parallel::mclapply(
    list_AitDist,function(x){
      choosing_alpha_Smetric(
        D = x,
        list_normalized_dist = list_normalized_geo_abiotics_dists,
        list_dist = list_geo_abiotics_dists,
        grid_alpha=seq(0,1,0.01),
        K=5,
      )
    },mc.cores = 10))

saveRDS(df_evaluation,paste0(savingdir,'/','df_evaluation_choosing_alpha'))



df_summary = bind_rows(
  df_evaluation %>% 
    group_by(method,alpha,DistMetric) %>% 
    summarise_all(.funs = mean) %>% mutate(summary_metric='Mean') %>% ungroup() %>% 
    pivot_longer(cols = -c(method,alpha,DistMetric,summary_metric),names_to = 'Metric'),

  df_evaluation %>% 
    group_by(method,alpha,DistMetric) %>% 
    summarise_all(.funs = sd) %>% mutate(summary_metric='SD') %>% ungroup() %>% 
    pivot_longer(cols = -c(method,alpha,DistMetric,summary_metric),names_to = 'Metric')
)


library(ggplot2)
df_summary %>%
  pivot_wider(id_cols = method:Metric,names_from = summary_metric ) %>% 
  select(Metric) %>% distinct() %>% pull


df_summary %>%
  pivot_wider(id_cols = method:Metric,names_from = summary_metric ) %>%
  filter(Metric=='avg_within_between_dist2medoid') %>% 
  ggplot(aes(x=alpha,y=Mean,color=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_wrap(Metric~DistMetric,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')

df_summary %>%
  pivot_wider(id_cols = method:Metric,names_from = summary_metric ) %>%
  filter(Metric=='avg_dist_within_dist2medoid') %>% 
  ggplot(aes(x=alpha,y=Mean,color=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_wrap(Metric~DistMetric,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')

df_summary %>%
  pivot_wider(id_cols = method:Metric,names_from = summary_metric ) %>%
  filter(Metric=='avg_dist_between_medoids') %>% 
  ggplot(aes(x=alpha,y=Mean,color=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_wrap(Metric~DistMetric,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')


df_summary %>%
  pivot_wider(id_cols = method:Metric,names_from = summary_metric ) %>%
  filter(Metric=='max_max_dist_within') %>% 
  ggplot(aes(x=alpha,y=Mean,color=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_wrap(Metric~DistMetric,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')

df_summary %>%
  pivot_wider(id_cols = method:Metric,names_from = summary_metric ) %>%
  filter(Metric=='avg_max_dist_within') %>% 
  ggplot(aes(x=alpha,y=Mean,color=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_wrap(Metric~DistMetric,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')


##### K= 8
df_evaluation <- data.table::rbindlist(
  parallel::mclapply(
    list_AitDist,function(x){
      choosing_alpha_Smetric(
        D = x,
        list_normalized_dist = list_normalized_geo_abiotics_dists,
        list_dist = list_geo_abiotics_dists,
        grid_alpha=seq(0,1,0.01),
        K=8,
      )
    },mc.cores = 10))

saveRDS(df_evaluation,paste0(savingdir,'/','df_evaluation_choosing_alpha_alt'))



df_summary = bind_rows(
  df_evaluation %>% 
    group_by(method,alpha,DistMetric) %>% 
    summarise_all(.funs = mean) %>% mutate(summary_metric='Mean') %>% ungroup() %>% 
    pivot_longer(cols = -c(method,alpha,DistMetric,summary_metric),names_to = 'Metric'),
  
  df_evaluation %>% 
    group_by(method,alpha,DistMetric) %>% 
    summarise_all(.funs = sd) %>% mutate(summary_metric='SD') %>% ungroup() %>% 
    pivot_longer(cols = -c(method,alpha,DistMetric,summary_metric),names_to = 'Metric')
)


library(ggplot2)
df_summary %>%
  pivot_wider(id_cols = method:Metric,names_from = summary_metric ) %>% 
  select(Metric) %>% distinct() %>% pull


df_summary %>%
  pivot_wider(id_cols = method:Metric,names_from = summary_metric ) %>%
  filter(Metric=='avg_within_between_dist2medoid') %>% 
  ggplot(aes(x=alpha,y=Mean,color=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_wrap(Metric~DistMetric,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')

df_summary %>%
  pivot_wider(id_cols = method:Metric,names_from = summary_metric ) %>%
  filter(Metric=='avg_dist_within_dist2medoid') %>% 
  ggplot(aes(x=alpha,y=Mean,color=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_wrap(Metric~DistMetric,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')

df_summary %>%
  pivot_wider(id_cols = method:Metric,names_from = summary_metric ) %>%
  filter(Metric=='avg_dist_between_medoids') %>% 
  ggplot(aes(x=alpha,y=Mean,color=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_wrap(Metric~DistMetric,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')


df_summary %>%
  pivot_wider(id_cols = method:Metric,names_from = summary_metric ) %>%
  filter(Metric=='max_max_dist_within') %>% 
  ggplot(aes(x=alpha,y=Mean,color=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_wrap(Metric~DistMetric,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')

df_summary %>%
  pivot_wider(id_cols = method:Metric,names_from = summary_metric ) %>%
  filter(Metric=='avg_max_dist_within') %>% 
  ggplot(aes(x=alpha,y=Mean,color=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.25)+
  facet_wrap(Metric~DistMetric,scales = 'free')+
  theme_minimal()+
  theme(legend.position = 'bottom')


