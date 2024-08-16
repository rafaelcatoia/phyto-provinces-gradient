########################################################################################################################################
## 01 - Generating Files:
## - data frame containing lat depth abiotics and long hurst provinces for each sample
## - list of distance matrices (Biotic, Abiotic, and Geographical)
## - list of distance matrices normalized (Biotic, Abiotic, and Geographical)
## - list containing perturbed versions of the Biotic component - to find k (number of clusters) and alpha (geographical structure)
########################################################################################################################################

### packages  ----------------------------------------------------------------------
library(dplyr) ; library(ggplot2)

## functions  ---------------------------------------------------------------------- 
### to filter ASVs
`%not_in%` <- purrr::negate(`%in%`)

### normalizing matrices
normalizeMatrix <- function(XX){
  normMat = norm(XX,type='2')
  return(XX/normMat)
}


### Getting data & path
root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
savingdir = root$find_file("saved_files")
path_df = root$find_file('data/phyto_gradients.csv') # filtered version

## filtering only the cruises that we want and depth
grump_longer <- data.table::fread(path_df)

grump_longer <- grump_longer %>%
  mutate(Raw.Sequence.Counts=Corrected_sequence_counts,
         ID_ASV=ASV_hash,
         SampleID=factor(SampleID))

### -----------------------------------------------------------
## Creating abiotics dataframe  -------------------------------
### -----------------------------------------------------------

vet_abiotic = c("Temperature","Salinity","Oxygen",
  "Silicate","NO2","NO3","PO4"
)

## --- we must correct the missing on the abiotics !
df_geo_abiotics <- grump_longer %>%
  select(SampleID,one_of(vet_abiotic),
         Latitude,Longitude,Depth,Longhurst_Short) %>%
  distinct() %>% arrange(SampleID)

saveRDS(df_geo_abiotics,file = paste0(savingdir,'/','df_geo_abiotics'))
### --------------------------------------------------------------
### Creating the list containing distance matrices 
### Normalized an unnormalized
### --------------------------------------------------------------
  

###### First the normalized ones  ################################
geoDist = df_geo_abiotics %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  transmute(lat_scaled = (Latitude - mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix() #%>% normalizeMatrix()

abioticDist = df_geo_abiotics %>% 
  transmute(
    Temperature = scale(Temperature),
    Salinity = scale(Salinity),
    Oxygen = scale(Oxygen),
    #Silicate = scale(Silicate),
    #NO2 = scale(NO2),
    #NO2 = scale(NO3),
    #PO4 = scale(PO4)
    ) %>% 
  as.matrix() %>% dist() %>% as.matrix() #%>%  normalizeMatrix()

min_raw_count=0.001

bioticDist = grump_longer %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,ID_ASV,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,ID_ASV) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = ID_ASV ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()


list_abio_bio_geo_dist <- list(
  geoDist = geoDist,
  abioticDist = abioticDist,
  bioticDist = bioticDist
  )

list_abio_bio_geo_dist_normallized <- list(
  geoDist = geoDist %>% normalizeMatrix(),
  abioticDist = abioticDist %>% normalizeMatrix(),
  bioticDist = bioticDist %>% normalizeMatrix()
)

##list_abio_bio_geo_dist contains all distance matrices
saveRDS(list_abio_bio_geo_dist,file = paste0(savingdir,'/','list_abio_bio_geo_dist'))

##list_abio_bio_geo_dist_normallized contains all distance matrices but they are normalized
saveRDS(list_abio_bio_geo_dist_normallized,file = paste0(savingdir,'/','list_abio_bio_geo_dist_normallized'))


### ----------------------------------------------------------------------
### Creating the B replicates of aitchison distance -- To tune K and alpha
### ----------------------------------------------------------------------

## number of replicates 
B=500

## list with perturbed versions of the original data
list_AitDist = list()

## Creating linst of all ASVs
idASVs = grump_longer %>% select(ID_ASV) %>% distinct() %>% pull

## percentage of columns
pct_colSubSample = 0.5

## minimal value to add in order to get the CLR transform / Aitchison distance
min_raw_count = grump_longer %>% select(Raw.Sequence.Counts) %>% min()
min_raw_count = min_raw_count/100

## In this case I'm fixing 0.001
min_raw_count = 0.001


## initiating the counter we will use a while loop, because when removing some ASVs
## it is possible that some samples are going to be excluded, and we don't
## want that to happen.

accepted_sample <- 1 
seedI = 5757

## untill we have less then or equal to B replicates:
while(accepted_sample<=B){
  ## sampling the ASVs to be used
  asvSubset=sample(idASVs,size = pct_colSubSample*length(idASVs))
  
  ## creating the aitchison distance matrix
  mat_repB <- grump_longer %>%
    ## first, filter the asvs that are in sample for this iteration
    filter(ID_ASV %in% asvSubset) %>% 
    
    ## adding a small "count" to every ASVs
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    
    ## now groupping by sample ID, ID_ASV
    group_by(SampleID,ID_ASV) %>% 
    
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    
    ## pivoting the table so we can create a matrix,
    ## where samples are in the rows and ASVs are in the columns
    tidyr::pivot_wider(id_cols = SampleID,names_from = ID_ASV ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    ## changing to data frame to avoid compatibility issues with other funcitos
    data.frame() %>% 
    
    ## here we make the rows to add to 1. 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    
    ## making sure sample ID is the first column
    relocate(SampleID,sort(names(.))) %>% 
    
    ## arranging so this distance matrix is in the same order
    arrange(SampleID) %>% 
    
    ## remove column sample id so we can calculate aitchison distance
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% 
    
    ## normalizing matrix
    normalizeMatrix()
  
  ## we repeate all the steps for the asvs that are not in the "in sample"
  
  mat_evalB <- grump_longer %>%
    filter(ID_ASV %not_in% asvSubset) %>% 
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  
    group_by(SampleID,ID_ASV) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = ID_ASV ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix()
  
  
  ## Here we check if both, the "In sample" matrix and the "out of sample" matrix have the right dimension
  if( sum(c(dim(mat_repB),dim(mat_evalB))!=nrow(df_geo_abiotics))  == 0){
    ## If both have the right dimension, we append the list.
    list_AitDist = rlist::list.append(list_AitDist,list(INS=mat_repB,OOS=mat_evalB))
    ## Updating the counter
    accepted_sample = accepted_sample+1
  } # if it is not accepted, than we move to the next iteration
  cat(paste0('Iteration number ------ ::',accepted_sample,'\n'))
}

## this list is used to calculate the metrics for choosing k and alpha.

saveRDS(list_AitDist,file = paste0(savingdir,'/','list_AitDist_IS_OOS'))

### -------------------------------------------------------------------------
### Creating the B replicates of Aitchison distance -- To find Provinces
### -------------------------------------------------------------------------
## number of replicates 
B=100

## list with perturbed versions of the original data
list_AitDist = list()

## Creating linst of all ASVs
idASVs = grump_longer %>% select(ID_ASV) %>% distinct() %>% pull

## percentage of columns
pct_colSubSample = 0.75

## minimal value to add in order to get the CLR transform / Aitchison distance
min_raw_count = grump_longer %>% select(Raw.Sequence.Counts) %>% min()
min_raw_count = min_raw_count/100
min_raw_count = 0.001
## initiating the counter we will use a while loop, because when removing some ASVs
## it is possible that some samples are going to be excluded, and we don't
## want that to happen.

accepted_sample <- 1 
seedI = 5757

while(accepted_sample<=B){
  ## sampling the ASVs to be used
  asvSubset=sample(idASVs,size = pct_colSubSample*length(idASVs))
  ## creating the aitchison distance matrix 
  mat_repB <- grump_longer %>%
    ## 
    filter(ID_ASV %in% asvSubset) %>% 
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    select(SampleID,ID_ASV,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,ID_ASV) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = ID_ASV ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% 
    normalizeMatrix()
  
  if( sum(dim(mat_repB) != nrow(df_geo_abiotics))  == 0){
    list_AitDist = rlist::list.append(list_AitDist,mat_repB)
    accepted_sample = accepted_sample+1
  }
  cat(paste0('Iteration number ------ ::',accepted_sample,'\n'))
}

saveRDS(list_AitDist,file = paste0(savingdir,'/','list_AitDist_IS_provinces'))

### -------------------------------------------------------------------------
### Now, moove to script 02! 
### -------------------------------------------------------------------------
