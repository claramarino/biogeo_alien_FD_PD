# Compute functional space for all bird species on 7 traits 
# and functional richness on its three first axes

rm(list = ls())

library(mFD)
library(ape)
library(dplyr)
library(picante)

func_traits <- readRDS("Data/03_func_traits.rds")

##### Functional space #####
traits <- func_traits[,c(1,3,5:10)]
rownames(traits) <- func_traits$Species1
traits$Trophic.Level<- factor(traits$Trophic.Level,order=TRUE,levels=c("Herbivore","Omnivore","Carnivore", "Scavenger"))
traits$Primary.Lifestyle<- as.factor(traits$Primary.Lifestyle)
traits$Mass <- log(traits$Mass)
traits$Beak.Length_Nares <- log(traits$Beak.Length_Nares)
traits$Beak.Depth <- log(traits$Beak.Depth)
traits <- na.omit(traits)

traits_cat <- data.frame( matrix(c("Beak.Length_Nares", "Q",
                                   "Beak.Depth", "Q", 
                                   "Hand.Wing.Index", "Q", 
                                   "Mass", "Q", 
                                   "Trophic.Level", "O", 
                                   "Primary.Lifestyle", "N", 
                                   "habitat.breadth", "Q"), 
                                 nrow = 7, ncol=2, byrow = T, 
                                 dimnames = list(rep("total",7),c("trait_name", "trait_type"))))


# Functional distance based on the traits
sp_dist <- funct.dist(sp_tr         = traits[,2:8],
                      tr_cat        = traits_cat,
                      metric        = "gower",
                      scale_euclid  = "scale_center",
                      ordinal_var   = "classic",
                      weight_type   = "equal",
                      stop_if_NA    = T)

saveRDS(sp_dist, "Data/13_sp_dist_world.rds")

# PCoA on the distance matrix
pcoa_out <- pcoa(sp_dist)
mat_coord<-pc_axes$vectors[,1:10]

eig <- pcoa_out$values$Eigenvalues/sum(pcoa_out$values$Eigenvalues[pcoa_out$values$Eigenvalues>0]) 
saveRDS(pcoa_out, "Data/13_pcoa_output_list_ordo.rds")
saveRDS(mat_coord, "Data/13_coord_10_FD_ordo.rds")

# Distances recalculated from the 3 first axes of the PCoA
traits_cat <- data.frame( matrix(c("Axis.1", "Q", 
                                   "Axis.2", "Q", 
                                   "Axis.3", "Q"), nrow = 3, ncol=2, byrow = T, dimnames = list(rep("total",3), c("trait_name", "trait_type"))))

new <- as.data.frame(mat_coord[,1:3])
sp_dist_pcoa <- funct.dist(sp_tr         = new,
                           tr_cat        = traits_cat,
                           metric        = "euclidean",
                           scale_euclid  = "scale_center",
                           ordinal_var   = "classic",
                           weight_type   = "equal",
                           stop_if_NA    = T)

dist <- dist(new)
Metrics::rmse(sp_dist_pcoa, dist)
# RMSE = 0.06291202

##### Functional richness (FRIC) ####

studied_isl <- read.csv("Data/Marino_et_al_DATA_407_ISL.csv", sep = ";")
site <-  readRDS("Data/01_species_isl_over3_matrix_name_traits_cor.rds")
site <- site[rownames(site)%in% studied_isl$ID, ]
pcoa_out <- readRDS("Data/13_pcoa_output_list_ordo.rds")
mat_coord <- readRDS("Data/13_coord_10_FD_ordo.rds")

# Changing names to fit with AVONET 
species_names <- unique(read.csv("Raw_Data/GAVIA_main_data_table.csv")$Binomial) 
species_names <- gsub(" ", "_", species_names)
Raw_AVONET <- read.csv("Raw_Data/AVONET_Raw_Data.csv")
Raw_AVONET$Species3_BirdTree <- sub(" ", "_", Raw_AVONET$Species3_BirdTree)
Raw_AVONET$Species1_BirdLife <- sub(" ", "_", Raw_AVONET$Species1_BirdLife)

tt <- func_traits[func_traits$Species1 %in% species_names,]
`%notin%` <- Negate(`%in%`)
pb <- species_names[species_names%notin% tt$Species1] 
to_change <- unique(Raw_AVONET[Raw_AVONET$Species3_BirdTree  %in% pb,c("Species1_BirdLife", "Species3_BirdTree")])
to_change <- to_change[to_change$Species1_BirdLife %in% func_traits$Species1,]
x <- species_names
for (i in (1:length(x))){
  if (x[i]%in%to_change$Species3_BirdTree){
    x[i] <- to_change[to_change$Species3_BirdTree==x[i],"Species1_BirdLife"]
  }
}

species_names <- x # 18 species not found in the trait database
traits_alien <- traits[traits$Species1 %in% species_names, ]

# Formatting community matrix
tot <- rep(1, dim(site)[2])
site <- rbind(site, tot)
comm <- matrix(c(1), nrow = 1, ncol = dim(traits_alien)[1])
colnames(comm)<- rownames(traits_alien)
comm <- plyr::rbind.fill(as.data.frame(site),  as.data.frame(comm))
comm[is.na(comm)] <- 0
rownames(comm)[1:407]<-rownames(site)
comm2 <- as.matrix(comm[1:408,])

# Functional diversity calculated on the three first axes of the PCoA
alpha_fd_indices <- alpha.fd.multidim(
  sp_faxes_coord   = mat_coord[,c("Axis.1", "Axis.2", "Axis.3")],
  asb_sp_w         = comm2[1:407,],
  ind_vect         = c("fric"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

saveRDS(alpha_fd_indices, "Data/13_func_div_world.rds")


# Null models for SES FD using all ALIEN birds 
null_fric <- as.data.frame(matrix(nrow = 407, ncol = 100))
for(j in 1:100){
  null_mod <- randomizeMatrix(comm2[1:407,], null.model = "richness")
  print(j)
  fric <- alpha.fd.multidim(
    sp_faxes_coord   = mat_coord[,c("Axis.1", "Axis.2", "Axis.3")],
    asb_sp_w         = null_mod,
    ind_vect         = c("fric"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  null_fric[,j] <- fric$functional_diversity_indices$fric 
  saveRDS(null_fric, "Data/13_null_model_fric_world_alien_only.rds")
}

func_div <- readRDS("Data/13_func_div_world.rds")

world_f_div <- data.frame(isl=rownames(site[1:407,]), 
                          SR=func_div$functional_diversity_indices$sp_richn,
                          fric=func_div$functional_diversity_indices$fric, 
                          SES_fric=(func_div$functional_diversity_indices$fric-
                                      rowMeans(null_fric))/matrixStats::rowSds(as.matrix(null_fric)))

saveRDS(world_f_div, "Data/13_func_div_ses_world_alien_only.rds")


# Null models for SES FD using all birds
null_fric <- as.data.frame(matrix(nrow = 1620, ncol = 100))
for(j in 1:100){
  null_mod <- randomizeMatrix(comm2[1:1620,], null.model = "richness")
  fric <- alpha.fd.multidim(
    sp_faxes_coord   = mat_coord[,c("Axis.1", "Axis.2", "Axis.3")],
    asb_sp_w         = null_mod,
    ind_vect         = c("fric"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  null_fric[,j] <- fric$functional_diversity_indices$fric 
  saveRDS(null_fric, "Data/13_null_model_fric_world.rds")
}

world_f_div <- data.frame(isl=rownames(site[1:1620,]), 
                          SR=func_div$functional_diversity_indices$sp_richn,
                          fric=func_div$functional_diversity_indices$fric, 
                          SES_fric=(func_div$functional_diversity_indices$fric-
                                      rowMeans(null_fric))/matrixStats::rowSds(as.matrix(null_fric)))

saveRDS(world_f_div, "Data/13_func_div_ses_world.rds")


##### Functional richness Natives #####
site <-  readRDS("Data/05_natives_isl_over3_matrix.rds")
site <- site[,colnames(site) %in% rownames(mat_coord)]

tot <- rep(1, dim(site)[2])
site <- rbind(site, tot)
comm <- matrix(c(1), nrow = 1, ncol = dim(traits)[1])
colnames(comm)<- rownames(traits)
comm <- plyr::rbind.fill(as.data.frame(site),  as.data.frame(comm))
comm[is.na(comm)] <- 0
rownames(comm)[1:1621]<-rownames(site)
comm2 <- as.matrix(comm[1:1621,])

pcoa_out <- readRDS("Data/13_pcoa_output_list_ordo.rds")
mat_coord <- readRDS("Data/13_coord_10_FD_ordo.rds")

Sys.time()
alpha_fd_indices <- alpha.fd.multidim(
  sp_faxes_coord   = mat_coord[,c("Axis.1", "Axis.2", "Axis.3")],
  asb_sp_w         = comm2[1:1620,],
  ind_vect         = c("fric"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)
Sys.time()

saveRDS(alpha_fd_indices, "Output/13_fric_natives_isl_over3.rds")

null_fric <- as.data.frame(matrix(nrow = 1620, ncol = 100))

for(j in 1:100){
  null_mod <- randomizeMatrix(comm2[1:1620,], null.model = "richness")
  fric <- alpha.fd.multidim(
    sp_faxes_coord   = mat_coord[,c("Axis.1", "Axis.2", "Axis.3")],
    asb_sp_w         = null_mod,
    ind_vect         = c("fric"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  null_fric[,j] <- fric$functional_diversity_indices$fric 
  saveRDS(null_fric, "Data/13_null_model_fric_natives.rds")
}

null_fric <- readRDS("Data/13_null_model_fric_natives.rds")
func_div <- readRDS("Output/13_fric_natives_isl_over3.rds")

nat_f_div <- data.frame(isl=rownames(site[1:1620,]), 
                        fric_nat=func_div$functional_diversity_indices$fric,
                        SES_fric_nat=(func_div$functional_diversity_indices$fric-
                                        rowMeans(null_fric))/matrixStats::rowSds(as.matrix(null_fric)))

saveRDS(nat_f_div, "Output/13_fric_+_ses_natives_isl_over3.rds")



#### Correlation traits and functional axes -------

pcoa_out <- readRDS("Data/13_pcoa_output_list_ordo.rds")
mat_coord <- readRDS("Data/13_coord_10_FD_ordo.rds")

plot(pcoa_out$vectors[,2], pc_axes$vectors[,1], type = "n", xlab = "PCoA1", ylab = "PCoA2",
     axes = TRUE, main = "PCoA sur les 10858 espèces")
text(pcoa_out$vectors[,2], pc_axes$vectors[,1], labels(func_traits$Species1), 
     cex = 0.9, xpd = TRUE)

eig <- pcoa_out$values$Eigenvalues/sum(pcoa_out$values$Eigenvalues[pcoa_out$values$Eigenvalues>0]) 
sum(eig[1:3])
eig[1:10]
pcoa_out$values


cor(pcoa_out$vectors[,1], traits$Beak.Length_Nares)

cor(pcoa_out$vectors[,1], as.numeric(traits$Trophic.Level), method = "spearman")
cor(pcoa_out$vectors[,2], as.numeric(traits$Trophic.Level), method = "spearman")
cor(pcoa_out$vectors[,3], as.numeric(traits$Trophic.Level), method = "spearman")

traits$Aerial <- 0
traits[traits$Primary.Lifestyle=="Aerial", ]$Aerial <- 1
cor(pcoa_out$vectors[,1], traits$Aerial, method = "spearman")
cor(pcoa_out$vectors[,2], traits$Aerial, method = "spearman")
cor(pcoa_out$vectors[,3], traits$Aerial, method = "spearman")

traits$Aquatic <- 0
traits[traits$Primary.Lifestyle=="Aquatic", ]$Aquatic <- 1
cor(pcoa_out$vectors[,1], traits$Aquatic, method = "spearman")
cor(pcoa_out$vectors[,2], traits$Aquatic, method = "spearman")
cor(pcoa_out$vectors[,3], traits$Aquatic, method = "spearman")

traits$Generalist <- 0
traits[traits$Primary.Lifestyle=="Generalist", ]$Generalist <- 1
cor(pcoa_out$vectors[,1], traits$Generalist, method = "spearman")
cor(pcoa_out$vectors[,2], traits$Generalist, method = "spearman")
cor(pcoa_out$vectors[,3], traits$Generalist, method = "spearman")

traits$Insessorial <- 0
traits[traits$Primary.Lifestyle=="Insessorial", ]$Insessorial <- 1
cor(pcoa_out$vectors[,1], traits$Insessorial, method = "spearman")
cor(pcoa_out$vectors[,2], traits$Insessorial, method = "spearman")
cor(pcoa_out$vectors[,3], traits$Insessorial, method = "spearman")

traits$Terrestrial <- 0
traits[traits$Primary.Lifestyle=="Terrestrial", ]$Terrestrial <- 1
cor(pcoa_out$vectors[,1], traits$Terrestrial, method = "spearman")
cor(pcoa_out$vectors[,2], traits$Terrestrial, method = "spearman")
cor(pcoa_out$vectors[,3], traits$Terrestrial, method = "spearman")

