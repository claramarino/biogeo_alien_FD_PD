# Compute SEM Models for 407 world oceanic islands 
# with at least 4 alien bird species

rm(list=ls())

library(lavaan)
library(piecewiseSEM)
library(lme4)
library(semEff)
library(tidyverse)

################ Load data ###############

# Clean data were obtained with script 00_create_clean_database
# Load final dataset with all scaled variables used in SEMs 
all_ok <- read.csv2("Data/Marino_et_al_DATA_407_ISL.csv")


# For a subset of islands (n=96), we could obtain data regarding
# the colonization pressure (CP), that is the total number of species
# introduced to an area (Blackburn et al., 2020)

# We thus perform SEM on those 96 islands with and without the CP variable 
# to check how this factor can influence our results.

# Filter the islands with CP data
isl_CP <- all_ok %>% filter(!is.na(CP))
table(isl_CP$Realm)
table(all_ok$Realm)


###### First models: without CP in the subset of 96 islands ######

pop <- lm(pop ~ Area + Lat + SLMP + Dist , data = isl_CP)
GDP <- lm(GDP ~ Area + pop + Dist + Elev + Lat, data = isl_CP)
connect <- lm(connect ~ Area + Lat + SLMP + Dist + pop, data = isl_CP)
static <- lm(static_modif ~ GDP + pop + modif_change + Area + Lat + connect, data = isl_CP)
modif <- lm(modif_change ~ GDP + pop + Elev + SLMP + connect, data = isl_CP)
# CP <- lm(CP ~ pop + GDP + connect + Dist + Area + static_modif, data = isl_CP)

nat_SR <- lm(SR_nat ~ Area + Dist + Lat + Elev + SLMP, data = isl_CP)
nat_PD <- lm(PD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat, data = isl_CP)
nat_FD <- lm(FD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat, data = isl_CP)

alien_SR <- lm(SR_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                 pop + GDP + modif_change + static_modif + connect, data = isl_CP)
alien_PD <- lm(PD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + PD_nat +
                 pop + GDP + modif_change + static_modif + connect + SR_alien, 
               data = isl_CP)
alien_FD <- lm(FD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + FD_nat +
                   pop + GDP + modif_change + static_modif + connect + SR_alien, 
                 data = isl_CP)

# build total model
# one for FD and one for pd

pd_mod <- psem(pop, GDP, connect, static, modif,
                  nat_SR, nat_PD, alien_SR, alien_PD,
                  SR_nat %~~% connect,
                  #SR_nat %~~% CP,
                  PD_nat %~~% pop,
                  PD_nat %~~% connect,
                  PD_nat %~~% GDP,
                  data = isl_CP)

fisherC(pd_mod) # F.C = 43.6, df = 34, P = 0.12 
rsquared(pd_mod)

FD_mod <- psem(pop, GDP, connect, static, modif, 
                    nat_SR, nat_FD, alien_SR, alien_FD,
                    SR_nat %~~% connect,
                    FD_nat %~~% pop,
                    FD_nat %~~% static_modif, 
                    FD_nat %~~% connect,
                    data = isl_CP)
d1 <- dSep(FD_mod)
fisherC(FD_mod) # F.C = 36, df = 34, P = .345 
rsquared(FD_mod)


###### Second models: with CP in the SEM ######

pop <- lm(pop ~ Area + Lat + SLMP + Dist , data = isl_CP)
GDP <- lm(GDP ~ Area + pop + Dist + Elev + Lat, data = isl_CP)
connect <- lm(connect ~ Area + Lat + SLMP + Dist + pop, data = isl_CP)
static <- lm(static_modif ~ GDP + pop + modif_change + Area + Lat + connect, data = isl_CP)
modif <- lm(modif_change ~ GDP + pop + Elev + SLMP + connect, data = isl_CP)
CP <- lm(CP ~ pop + GDP + connect + Dist + Area + static_modif, data = isl_CP)
# biotic context depends on biogeographic variables
# native SR drives native FD & PD
nat_SR <- lm(SR_nat ~ Area + Dist + Lat + Elev + SLMP, data = isl_CP)
nat_PD <- lm(PD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat, data = isl_CP)
nat_FD <- lm(FD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat, data = isl_CP)
# alien diversities depend on all contextual variables
alien_SR <- lm(SR_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                   pop + GDP + modif_change + static_modif + connect + CP, data = isl_CP)
alien_PD <- lm(PD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + PD_nat +
                   pop + GDP + modif_change + static_modif + connect + CP + SR_alien, 
                 data = isl_CP)
alien_FD <- lm(FD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + FD_nat +
                     pop + GDP + modif_change + static_modif + connect + CP + SR_alien, 
                   data = isl_CP)

# build total model
# one for FD and one for pd

cp_pd_mod <- psem(pop, GDP, connect, static, modif, CP,
               nat_SR, nat_PD, alien_SR, alien_PD,
               SR_nat %~~% connect,
               SR_nat %~~% CP,
               PD_nat %~~% pop,
               PD_nat %~~% connect,
               PD_nat %~~% GDP,
               data = isl_CP)

#basisSet(cp_pd_mod)
#d1 <- dSep(cp_pd_mod)
fisherC(cp_pd_mod) 
#d1[which.max(abs(d1$Crit.Value)),] 

cp_FD_mod <- psem(pop, GDP, connect, static, modif, CP,
                 nat_SR, nat_FD, alien_SR, alien_FD,
                 SR_nat %~~% connect,
                 SR_nat %~~% CP,
                 FD_nat %~~% pop,
                 FD_nat %~~% static_modif, 
                 FD_nat %~~% connect,
                 data = isl_CP)
#d1 <- dSep(cp_FD_mod)
fisherC(cp_FD_mod)
#d1[which.max(abs(d1$Crit.Value)),] 

rsquared(cp_FD_mod)
rsquared(cp_pd_mod)

########## Get direct/ indirect effects of the 2 models without CP  ######

R_sample = 10000

Sys.time()

boot_mod_fd <- bootEff(FD_mod, R = R_sample)
eff_fd <- semEff(boot_mod_fd)

Sys.time()

boot_mod_pd <- bootEff(pd_mod, R = R_sample)
eff_pd <- semEff(boot_mod_pd)

Sys.time()

##### Save outputs  ######

# fd
saveRDS(FD_mod, "Output/11_SEM_noCP_FD_model.rds")
saveRDS(eff_fd, "Output/11_SEM_noCP_FD_boot_eff.rds")
# pd
saveRDS(pd_mod, "Output/11_SEM_noCP_PD_model.rds")
saveRDS(eff_pd, "Output/11_SEM_noCP_PD_boot_eff.rds")


########## Get direct/ indirect effects of the 2 models with CP  ######

R_sample = 10000

Sys.time()

cp_boot_mod_fd <- bootEff(cp_FD_mod, R = R_sample)
cp_eff_fd <- semEff(cp_boot_mod_fd)

Sys.time()

cp_boot_mod_pd <- bootEff(cp_pd_mod, R = R_sample)
cp_eff_pd <- semEff(cp_boot_mod_pd)

Sys.time()

##### Save outputs 

# fd
saveRDS(cp_FD_mod, "Output/11_SEM_CP_FD_model.rds")
saveRDS(cp_eff_fd, "Output/11_SEM_CP_FD_boot_eff.rds")
# pd
saveRDS(cp_pd_mod, "Output/11_SEM_CP_PD_model.rds")
saveRDS(cp_eff_pd, "Output/11_SEM_CP_PD_boot_eff.rds")





############## Simple model outputs

# for PD 
rsquared(cp_pd_mod)
coefs <- coefs(cp_pd_mod)
coefs[coefs$P.Value<0.05,]


# for FD 
rsquared(cp_FD_mod)
coefs <- coefs(cp_FD_mod) 
coefs[coefs$P.Value<0.05,]