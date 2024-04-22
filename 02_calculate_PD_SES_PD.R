# Phylogeny selection, PD and SES PD calculation for alien birds 
# and native assemblages 

rm(list = ls())

library(dplyr)
library(picante)
library(tidyverse)
library(ape)

#### Phylogenetic trees for PD calculation ####

# Pick randomly 1000 trees from Jeltz 2012
filenames <- list.files("Raw_Data/Trees/", pattern="*.tre", full.names=TRUE)
all_trees <- lapply(filenames, read.tree)
final_list <- list()
for(i in (1:10)){
  for(j in (1:1000)){
    final_list <- append(final_list, all_trees[[i]][j])
  }
}
# Delete 2 first incomplete lists
final_list <- final_list[3:10001]
set.seed(47)
phyloG <- sample(final_list, 100)

saveRDS(phyloG,"Data/02_trees.rds")
phyloG <- readRDS("Data/02_trees.rds")

# Alien communities matrix
site <-  readRDS("Data/01_species_isl_over3_matrix.rds")

# Renaming species absent from the phylogenetic trees
`%notin%` <- Negate(`%in%`)
spe <- as.data.frame(phyloG[[1]]$tip.label)
sb <- colnames(site)[colnames(site) %notin% spe$`phyloG[[1]]$tip.label`]
site[site[,"Acridotheres_javanicus"]==1,][,"Acridotheres_grandis"] <- 1
site[site[,"Streptopelia_risoria"]==1,][,"Streptopelia_roseogrisea"] <- 1
colnames(site)[colnames(site)=="Porphyrio_mantelli"]<-"Porphyrio_hochstetteri"
colnames(site)[colnames(site)=="White_Cockatoo"]<-"Cacatua_alba"

site <- as.matrix(select(as.data.frame(site), -c("Acridotheres_javanicus",
                                                 "Established",
                                                 "Streptopelia_risoria")))


#### Phylogenetic Diversity (PD) #### 

# Global PD of alien birds 
species <- matrix(c(1),nrow=1,ncol=dim(site)[2],byrow = FALSE)
colnames(species) <- colnames(site)

pd_glob <- data.frame()
for (i in (1:length(phyloG))){
  pd_glob <- rbind(pd_glob, pd(species, phyloG[[i]]))
  if (i %in% seq(0, 100, by=10)){
    saveRDS(pd_glob, "Data/11_pd_glob_all_isl.rds")
  }
}

# Null models for global PD based on all birds
species_names <- phyloG[[1]]$tip.label

mat <- matrix(c(1), nrow = 1, ncol = dim(site)[2], byrow = F)
null_glob_pd <- data.frame()
for (j in  (1:100)){
  names <- sample(species_names,dim(site)[2])
  for (i in (1:length(phyloG))){
    colnames(mat)<- names
    null_glob_pd <- rbind(null_glob_pd, pd(mat, phyloG[[i]]))
    print(i)
    if (i %in% seq(0, 1000, by=50)){
      saveRDS(null_glob_pd, "Data/11_null_pd_glob_all_isl.rds")
    }
  }
}

null_glob_pd <- readRDS("Data/11_null_pd_glob_all_isl.rds")

# Null models for global PD based on all ALIEN birds
species_names <- unique(read.csv("Raw_Data/GAVIA_main_data_table.csv")$Binomial) # World alien species from GAVIA
species_names <- gsub(" ", "_", species_names)
phyloG_name <- phyloG[[1]]$tip.label
species_names <- species_names[species_names %in% phyloG_name] # Removing 9 alien species not in Jetz phylogeny

mat <- matrix(c(1), nrow = 1, ncol = dim(site)[2], byrow = F)
null_glob_pd <- data.frame()
for (j in  (1:100)){
  names <- sample(species_names,dim(site)[2])
  for (i in (1:length(phyloG))){
    colnames(mat)<- names
    null_glob_pd <- rbind(null_glob_pd, pd(mat, phyloG[[i]]))
    print(i)
    if (i %in% seq(0, 1000, by=50)){
      saveRDS(null_glob_pd, "Data/11_null_pd_glob_alien_isl.rds")
    }
  }
}

null_glob_pd <- readRDS("Data/11_null_pd_glob_alien_isl.rds")


# PD by island 
all_pd <- data.frame()
all_pd <- pd(site, phyloG[[1]])
for (i in (2:length(phyloG))){
  all_pd[,i] <- pd(site, phyloG[[i]])
  if (i %in% seq(0, 100, by=10)){
    saveRDS(all_pd, "Data/11_pd_world.rds")
  }
}

# Null models for SES PD by island using all birds
null_pd <- as.data.frame(matrix(nrow = 1620*3, ncol = 100))

species_names <- phyloG[[1]]$tip.label
null_mod <- site[1:1620,]
set.seed(42)
phylo3 <- sample(phyloG, 3)

for(j in 1:100){
  names <- sample(species_names,dim(null_mod)[2])
  colnames(null_mod) <- names
  null_mod <- randomizeMatrix(null_mod, null.model = "independentswap")
  for(i in 1:length(phylo3)){
    null_pd[(1620*(i-1)+1):(1620*i),j] <- pd(null_mod, phylo3[[i]])[,1]
  }
  if (j %in% seq(0, 100, by=10)){
    saveRDS(null_pd, "Data/11_null_model_world.rds")
  }
}

null_pd <- readRDS("Data/11_null_model_world.rds")
all_pd <- readRDS("Data/11_pd_world.rds")

tot_null_pd <- data.frame()
tot_null_pd[1:1620,1:100] <- null_pd[1:1620,1:100]
tot_null_pd[1:1620,101:200]<- null_pd[1621:3240,1:100]
tot_null_pd[1:1620,201:300]<- null_pd[3241:4860,1:100]
rownames(tot_null_pd) <- rownames(site)

pd_by_isl <- data.frame(isl=rownames(site), SR= rowSums(site), 
                        PD_mean=rowMeans(all_pd), 
                        ses_PD = (rowMeans(all_pd)-rowMeans(tot_null_pd))/matrixStats::rowSds(as.matrix(tot_null_pd)))
saveRDS(pd_by_isl, "Data/11_pd_by_isl_world.rds")

# Null models for SES PD by island using all ALIEN birds
null_pd <- as.data.frame(matrix(nrow = 1620*3, ncol = 100))

species_names <- unique(read.csv("Raw_Data/GAVIA_main_data_table.csv")$Binomial) 
species_names <- gsub(" ", "_", species_names)
phyloG_name <- phyloG[[1]]$tip.label
species_names <- species_names[species_names %in% phyloG_name]
null_mod <- site[1:1620,]
set.seed(42)
phylo3 <- sample(phyloG, 3)

for(j in 1:100){
  names <- sample(species_names,dim(null_mod)[2])
  colnames(null_mod) <- names
  null_mod <- randomizeMatrix(null_mod, null.model = "independentswap")
  for(i in 1:length(phylo3)){
    null_pd[(1620*(i-1)+1):(1620*i),j] <- pd(null_mod, phylo3[[i]])[,1]
  }
  if (j %in% seq(0, 100, by=10)){
    saveRDS(null_pd, "Data/11_null_model_world_alien.rds")
  }
}

null_pd <- readRDS("Data/11_null_model_world_alien.rds")
all_pd <- readRDS("Data/11_pd_world.rds")

tot_null_pd <- data.frame()
tot_null_pd[1:1620,1:100] <- null_pd[1:1620,1:100]
tot_null_pd[1:1620,101:200]<- null_pd[1621:3240,1:100]
tot_null_pd[1:1620,201:300]<- null_pd[3241:4860,1:100]
rownames(tot_null_pd) <- rownames(site)

pd_by_isl <- data.frame(isl=rownames(site), SR= rowSums(site), 
                        PD_mean=rowMeans(all_pd), 
                        ses_PD = (rowMeans(all_pd)-rowMeans(tot_null_pd))/matrixStats::rowSds(as.matrix(tot_null_pd)))

saveRDS(pd_by_isl, "Data/11_pd_by_isl_world_alien_only.rds")

##### Phylogenetic diversity natives ####
library(rgbif)
site <-  readRDS("Data/05_natives_isl_over3_matrix.rds")

# Renaming species not in the phylogenetic trees
`%notin%` <- Negate(`%in%`)
spe <- as.data.frame(phylogeny[[1]]$tip.label)
pb <- colnames(site)[colnames(site)%notin% spe$`phylogeny[[1]]$tip.label`] 
Raw_AVONET <- read.csv("Raw_Data/AVONET_Raw_Data.csv")
Raw_AVONET$Species3_BirdTree <- sub(" ", "_", Raw_AVONET$Species3_BirdTree)
Raw_AVONET$Species1_BirdLife <- sub(" ", "_", Raw_AVONET$Species1_BirdLife)

to_change <- unique(Raw_AVONET[Raw_AVONET$Species1_BirdLife  %in% pb,c("Species1_BirdLife", "Species3_BirdTree")])
to_change <- to_change[to_change$Species3_BirdTree %in% spe$`phylogeny[[1]]$tip.label`,]

x <- colnames(site)
for (i in (1:length(x))){
  if (x[i] %in% to_change$Species1_BirdLife ){
    x[i] <- to_change[to_change$Species1_BirdLife==x[i],"Species3_BirdTree"][1]
    
  }
}
colnames(site) <- x
# Only Dicaeum dayakorum left : native from Borneo and described in 2019

# PD all islands
all_pd <- data.frame()
for (i in (2:length(phyloG))){
  all_pd[,i] <- pd(site, phyloG[[i]])
  print(i)
  if (i %in% seq(0, 100, by=10)){
    print("saving")
    saveRDS(all_pd, "Output/11_pd_natives_world.rds")
  }
}

null_pd <- as.data.frame(matrix(nrow = 1620*3, ncol = 100))
species_names <- phyloG[[1]]$tip.label
null_mod <- site[1:1620,]

set.seed(42)
phylo3 <- sample(phyloG, 3)

for(j in 1:100){
  names <- sample(species_names,dim(null_mod)[2])
  colnames(null_mod) <- names
  null_mod <- randomizeMatrix(null_mod, null.model = "independentswap")
  for(i in 1:length(phylo3)){
    null_pd[(1620*(i-1)+1):(1620*i),j] <- pd(null_mod, phylo3[[i]])[,1]
  }
  if (j %in% seq(0, 100, by=10)){
    saveRDS(null_pd, "Output/11_null_model_natives_world.rds")
  }
}

null_pd <- readRDS("Output/11_null_model_natives_world.rds")
all_pd <- readRDS("Output/11_pd_natives_world.rds")

tot_null_pd <- data.frame()
tot_null_pd[1:1620,1:100] <- null_pd[1:1620,1:100]
tot_null_pd[1:1620,101:200]<- null_pd[1621:3240,1:100]
tot_null_pd[1:1620,201:300]<- null_pd[3241:4860,1:100]

rownames(tot_null_pd)<- rownames(site)

pd_by_isl <- data.frame(isl=rownames(site), SR= rowSums(site), 
                        PD_mean=rowMeans(all_pd), 
                        ses_PD = (rowMeans(all_pd)-rowMeans(tot_null_pd))/matrixStats::rowSds(as.matrix(tot_null_pd)))

saveRDS(pd_by_isl, "Output/11_pd_+_ses_natives_world.rds")
