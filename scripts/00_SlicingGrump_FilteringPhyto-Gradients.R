########################################################
## Script that generates all the files of the shiny app
########################################################

## Packages ########################################################
library(dplyr) ; library(tidyr)

## Directories #############################################################
root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")

## Loading grump data  #####################################################
path_grump = root$find_file('data/grump_asv_long.csv') 
df_grump_all <- data.table::fread(path_grump)

## Selecting the slice that we are interested in  ##########################
## For this script we are only looking at the Gradients 2 and 3 ############
grump_slice = df_grump_all %>%
  filter(Cruise %in% c('Gradients_2','Gradients_3')) %>% 
  filter(Latitude > 25)

## Looking at the sample level 
vet_abiotic = c(
  "Temperature","Salinity","Oxygen",
  "Silicate","NO2","NO3","PO4")

## Computting missing on the abiotic factors
grump_slice  %>% select(
  SampleID,
  any_of(vet_abiotic),
  Latitude,Longitude,Depth,Longhurst_Short) %>%
  distinct() %>% arrange(SampleID) %>% 
  select(one_of(vet_abiotic)) %>% 
  apply(.,2,function(x){sum(is.na(x))})


## And we also want to look at only the phytoplanktons
grump_slice[, c(2:11,54)][is.na(grump_slice[, c(2:11,54)])] <- "Not Applicable" 
cyano <- rbind(grump_slice[grump_slice$Eco_relevant_plank_groups=="Prochlorococcus",],
               grump_slice[grump_slice$Eco_relevant_plank_groups=="Synechococcus",],
               grump_slice[grepl("Richelia", grump_slice$Genus),],
               grump_slice[grepl("Trichodesmium", grump_slice$Genus),],
               grump_slice[grepl("UCYN-A", grump_slice$Genus),],
               grump_slice[grepl("Crocosphaera", grump_slice$Genus),],
               grump_slice[grepl("UCYN-C", grump_slice$Genus),])
## Subset out the chloroplast 16S ASVs from the main dataframe
chl.16s <- grump_slice[grump_slice$Sequence_Type=="Chloroplast_16S",]
## Remove the following ASVs from that subset
chl.16s <- chl.16s[!chl.16s$Supergroup=="Rhizaria",]
chl.16s <- chl.16s[!chl.16s$Supergroup=="Excavata",]
chl.16s <- chl.16s[!chl.16s$Supergroup=="Alveolata",]
chl.16s <- chl.16s[!chl.16s$Division=="Rhodophyta",]
## Final "phytoplankton" community we will use for the time being

phyto.all <- rbind(cyano, chl.16s)

## Saving All ASVs but only for P16SN
data.table::fwrite(phyto.all,paste0(datadir,'/phyto_gradients.csv'))

