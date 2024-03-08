
rm(list=ls())

library(mFD)
library(ape)
library(dplyr)
library(picante)

func_traits <- readRDS("Data/03_func_traits.rds")
site <-  readRDS("Data/01_species_sites_matrix.rds")

##### Functionnal space #####
traits <- func_traits[,c(1,3,5:10)]
rownames(traits) <- func_traits$Species1
traits$Trophic.Level<- factor(traits$Trophic.Level,order=TRUE,levels=c("Herbivore","Omnivore","Carnivore", "Scavenger"))
traits$Primary.Lifestyle<- as.factor(traits$Primary.Lifestyle)
traits$Mass <- log(traits$Mass)
traits$Beak.Length_Nares <- log(traits$Beak.Length_Nares)
traits$Beak.Depth <- log(traits$Beak.Depth)
traits <- na.omit(traits)

##### Standard Fric #####
traits_cat <- data.frame( matrix(c("Beak.Length_Nares", "Q",
                                   "Beak.Depth", "Q", 
                                   "Hand.Wing.Index", "Q", 
                                   "Mass", "Q", 
                                   "Trophic.Level", "O", 
                                   "Primary.Lifestyle", "N", 
                                   "habitat.breadth", "Q"), nrow = 7, ncol=2, byrow = T, dimnames = list(rep("total",7), 
                                                                                                         c("trait_name", "trait_type"))))


sp_dist <- funct.dist(sp_tr         = traits[,2:8],
                      tr_cat        = traits_cat,
                      metric        = "gower",
                      scale_euclid  = "scale_center",
                      ordinal_var   = "classic",
                      weight_type   = "equal",
                      stop_if_NA    = T)
saveRDS(sp_dist, "Data/13_sp_dist_world_240220.rds")

pc_axes <- pcoa(sp_dist)
mat_coord<-pc_axes$vectors[,1:10]

saveRDS(pc_axes, "Data/13_pcoa_output_list_ordo_240220.rds")
saveRDS(mat_coord, "Data/13_coord_10_FD_ordo_240220.rds")

pcoa_out <- readRDS("Data/13_pcoa_output_list_ordo_240220.rds")
mat_coord <- readRDS("Data/13_coord_10_FD_ordo_240220.rds")

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

# Distances recalculated from 3 first axes
traits_cat <- data.frame( matrix(c("Axis.1", "Q", 
                                   "Axis.2", "Q", 
                                   "Axis.3", "Q"), nrow = 3, ncol=2, byrow = T, dimnames = list(rep("total",3), + c("trait_name", "trait_type"))))
new <- as.data.frame(mat_coord[,1:3])
sp_dist_pcoa <- funct.dist(sp_tr         = new,
                           tr_cat        = traits_cat,
                           metric        = "euclidean",
                           scale_euclid  = "scale_center",
                           ordinal_var   = "classic",
                           weight_type   = "equal",
                           stop_if_NA    = T)


dist <- dist(new)
Metrics::rmse(sp_dist, dist)
# RMSE = 0.06291202


##### Fric, Feve, Fdiv of the world alien communities ####
site <-  readRDS("Data/01_species_isl_over3_matrix_name_traits_cor.rds")

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


alpha_fd_indices <- alpha.fd.multidim(
  sp_faxes_coord   = mat_coord[,c("Axis.1", "Axis.2", "Axis.3")],
  asb_sp_w         = comm2[1:1620,],
  ind_vect         = c("fric", "feve", "fdiv"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

saveRDS(alpha_fd_indices, "Data/13_func_div_world_240220.rds")

null_fric <- as.data.frame(matrix(nrow = 1620, ncol = 100))
null_feve <- as.data.frame(matrix(nrow = 1620, ncol = 100))
null_fdiv <- as.data.frame(matrix(nrow = 1620, ncol = 100))

for(j in 1:100){
  null_mod <- randomizeMatrix(comm2[1:1620,], null.model = "richness")
  print("modèle")
  print(j)
  fric <- alpha.fd.multidim(
    sp_faxes_coord   = mat_coord[,c("Axis.1", "Axis.2", "Axis.3")],
    asb_sp_w         = null_mod,
    ind_vect         = c("fric", "feve", "fdiv"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  null_fric[,j] <- fric$functional_diversity_indices$fric 
  null_feve[,j] <- fric$functional_diversity_indices$feve
  null_fdiv[,j] <- fric$functional_diversity_indices$fdiv
  
  print("saving")
  saveRDS(null_fric, "Data/13_null_model_fric_world.rds")
  saveRDS(null_feve, "Data/13_null_model_feve_world.rds")
  saveRDS(null_fdiv, "Data/13_null_model_fdiv_world.rds")
  
}

world_f_div <- data.frame(isl=rownames(site[1:1620,]), 
                          SR=func_div$functional_diversity_indices$sp_richn,
                          fric=func_div$functional_diversity_indices$fric, 
                          feve=func_div$functional_diversity_indices$feve, 
                          fdiv=func_div$functional_diversity_indices$fdiv,
                          SES_fric=(func_div$functional_diversity_indices$fric-
                                      rowMeans(null_fric))/matrixStats::rowSds(as.matrix(null_fric)),
                          SES_feve=(func_div$functional_diversity_indices$feve-
                                      rowMeans(null_feve))/matrixStats::rowSds(as.matrix(null_feve)),
                          SES_fdiv=(func_div$functional_diversity_indices$fdiv-
                                      rowMeans(null_fdiv))/matrixStats::rowSds(as.matrix(null_fdiv)))

saveRDS(world_f_div, "Data/13_func_div_ses_world.rds")
world_f_div <- readRDS( "Data/13_func_div_ses_world.rds")
boxplot(world_f_div$SES_fric, xlab="FRic", ylab="SES")
mean(world_f_div$SES_fric)


# check if fd indices are the good ones
fdtoday <- readRDS("Data/13_func_div_world_240220.rds")

world_f_div <- data.frame(isl=rownames(site[1:1620,]), 
                          SR_new=fdtoday$functional_diversity_indices$sp_richn,
                          fric_new=fdtoday$functional_diversity_indices$fric)


fric_alien <- readRDS("Data/13_func_div_ses_world.rds")


test <- left_join(world_f_div, fric_alien %>% select(isl, SR, fric))

plot(test$SR_new, test$SR)
plot(test$fric_new, test$fric)

sum(round(test$fric_new,5) == round(test$fric,5))
