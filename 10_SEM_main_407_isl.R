# Compute SEM Models for 407 world oceanic islands 
# with at least 4 alien bird species

rm(list=ls())

library(lavaan)
library(piecewiseSEM)
library(lme4)
library(semEff)
library(tidyverse)
# library(ggpubr) # useful for what?

################ Load data ###############

# Clean data were obtained with script 00_create_clean_database
# Load final dataset with all scaled variables used in SEMs 
all_ok <- read.csv2("Data/Marino_et_al_DATA_407_ISL.csv")
colSums(is.na(all_ok))

# remove CP and TFSI for main models because of NA
all_ok <- all_ok %>% select(-c(ARCHIP, ISLAND, CP, TSFI))

###################### Model specification #########################

# Specify submodels with realm as random effect
colnames(all_ok)


# human variables depend on biogeographic var
# and are related to each other
pop <- lmer(pop ~ Area + Lat + SLMP + Dist +
              (1|Realm), data = all_ok)
GDP <- lmer(GDP ~ Area + pop + Dist + Elev +
              (1|Realm), data = all_ok)
connect <- lmer(connect ~ Area + Lat + SLMP + Dist + pop + 
                  (1|Realm), data = all_ok)
static <- lmer(static_modif ~ GDP + pop + modif_change + Area + Lat + connect + 
                 (1|Realm), data = all_ok)
modif <- lmer(modif_change ~ GDP + pop + Elev  + SLMP+ connect +
                (1|Realm), data = all_ok)
# biotic context depends on biogeographic variables
# native SR drives native fric & PD
nat_SR <- lmer(SR_nat ~ Area + Dist + Lat + Elev + SLMP + 
                 (1|Realm), data = all_ok)
nat_PD <- lmer(PD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                 (1|Realm), data = all_ok)
nat_fric <- lmer(FD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat +
                   (1|Realm), data = all_ok)
# alien diversities depend on all contextual variables
alien_SR <- lmer(SR_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                   pop + GDP + modif_change + static_modif + connect +
                   (1|Realm), data = all_ok)
alien_PD <- lmer(PD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + PD_nat +
                   pop + GDP + modif_change + static_modif + connect + SR_alien + 
                   (1|Realm), data = all_ok)
alien_fric <- lmer(FD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + FD_nat +
                     pop + GDP + modif_change + static_modif + connect + SR_alien +
                     (1|Realm), data = all_ok)

# build total model
# one for fric and one for pd

pd_mod <- psem(pop, GDP, connect, static, modif, 
               nat_SR, nat_PD, alien_SR, 
               alien_PD,
               SR_nat %~~% connect,
               PD_nat %~~% pop,
               PD_nat %~~% connect,
               PD_nat %~~% GDP,
               data = all_ok)

basisSet(pd_mod)
d1 <- dSep(pd_mod)
fisherC(pd_mod) # F.C = 43.6, df = 36, P = 0.179
d1[which.max(abs(d1$Crit.Value)),] 
# adding correlated error between native diversities and human factors
# until the pval is above .05 => ok here

fric_mod <- psem(pop, GDP, connect, static, modif, 
                 nat_SR, nat_fric, alien_SR, 
                 alien_fric,
                 SR_nat %~~% connect,
                 FD_nat %~~% pop,
                 FD_nat %~~% static_modif, 
                 FD_nat %~~% connect,
                 data = all_ok)
d1 <- dSep(fric_mod)
fisherC(fric_mod) # F.C = 41.49, df = 36, P = 0.243 
d1[which.max(abs(d1$Crit.Value)),] 
# adding correlated error between native diversities and human factors
# until the pval is above .05 => ok here


############## Simple model outputs #############

# for PD 
rsquared(pd_mod) # Fort effet sur la richesse spé native et alien cf. R² marginal/conditional
coefs <- coefs(pd_mod)
coefs[coefs$P.Value<0.05,]


# for fric 
rsquared(fric_mod) # Mm chose qu'avant + effet assez important sur Fric (R²+0.2)
coefs <- coefs(fric_mod) 
coefs[coefs$P.Value<0.05,]


########### Response variables: SES-PD and SES-FD ###################

# Hypotheses are the same as above but we take the SES as response variables
# and not the diversity values per se

# human variables depend on biogeographic var
# and are related to each other
pop <- lmer(pop ~ Area + Lat + SLMP + Dist +
              (1|Realm), data = all_ok)
GDP <- lmer(GDP ~ Area + pop + Dist + Elev +
              (1|Realm), data = all_ok)
connect <- lmer(connect ~ Area + Lat + SLMP + Dist + pop + 
                  (1|Realm), data = all_ok)
static <- lmer(static_modif ~ GDP + pop + modif_change + Area + Lat + connect + 
                 (1|Realm), data = all_ok)
modif <- lmer(modif_change ~ GDP + pop + Elev  + SLMP+ connect +
                (1|Realm), data = all_ok)
# biotic context depends on biogeographic variables
# native SR drives native fric & PD
nat_SR <- lmer(SR_nat ~ Area + Dist + Lat + Elev + SLMP + 
                 (1|Realm), data = all_ok)
nat_SES_PD <- lmer(SES_PD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                     (1|Realm), data = all_ok)
nat_SES_fric <- lmer(SES_FD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat +
                       (1|Realm), data = all_ok)
# alien diversities depend on all contextual variables
alien_SR <- lmer(SR_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                   pop + GDP + modif_change + static_modif + connect +
                   (1|Realm), data = all_ok)
alien_SES_PD <- lmer(SES_PD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + SES_PD_nat +
                       pop + GDP + modif_change + static_modif + connect + SR_alien + 
                       (1|Realm), data = all_ok)
alien_SES_fric <- lmer(SES_FD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + SES_FD_nat +
                         pop + GDP + modif_change + static_modif + connect + SR_alien +
                         (1|Realm), data = all_ok)

# build total model
# one for SES_PD and one for SES_fric
# alien diversities depend on all contextual variables
alien_SR2 <- lmer(SR_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + SES_PD_nat +
                    pop + GDP + modif_change + static_modif + connect +
                    (1|Realm), data = all_ok)
SES_pd_mod <- psem(pop, GDP, connect, static, modif, 
                   nat_SR, nat_SES_PD, alien_SR2, alien_SES_PD,
                   SR_nat %~~% connect,
                   SES_PD_nat %~~% pop,
                   SES_PD_nat %~~% connect,
                   SES_PD_nat %~~% GDP,
                   data = all_ok)

d1 <- dSep(SES_pd_mod)
fisherC(SES_pd_mod) # F.C = 41.702, df = 34, P = 0.171 
d1[which.max(abs(d1$Crit.Value)),] 
# adding correlated error between native diversities and human factors
# until the pval is above .05
# + SES PD nat as a predictor of alien SR

SES_fric_mod <- psem(pop, GDP, connect, static, modif, 
                     nat_SR, nat_SES_fric, alien_SR, alien_SES_fric,
                     SR_nat %~~% connect,
                     SES_FD_nat %~~% pop,
                     data = all_ok)
d1 <- dSep(SES_fric_mod)
fisherC(SES_fric_mod) # F.C = 51.509, df = 40, P = 0.105 
d1[which.max(abs(d1$Crit.Value)),] 
# adding correlated error between native diversities and human factors
# until the pval is above .05 (less additions than for Fric)



############## Simple model outputs #############

# for PD 
rsquared(SES_pd_mod) # Fort effet sur la richesse spé native et alien cf. R² marginal/conditional
coefs <- coefs(SES_pd_mod)
coefs[coefs$P.Value<0.05,]


# for fric 
rsquared(SES_fric_mod) # Mm chose qu'avant + effet assez important sur Fric (R²+0.2)
coefs <- coefs(SES_fric_mod) 
coefs[coefs$P.Value<0.05,]


########## Get direct/ indirect effects of the 4 models #############

#install.packages("semEff")

# resample for confidence intervals
R_sample = 10000

Sys.time()

boot_mod_fd <- bootEff(fric_mod, R = R_sample, ran.eff = "Realm")
eff_fd <- semEff(boot_mod_fd)

Sys.time()

boot_mod_pd <- bootEff(pd_mod, R = R_sample, ran.eff = "Realm")
eff_pd <- semEff(boot_mod_pd)

Sys.time()

boot_mod_fd_ses <- bootEff(SES_fric_mod, R = R_sample, ran.eff = "Realm")
eff_fd_ses <- semEff(boot_mod_fd_ses)

Sys.time()

boot_mod_pd_ses <- bootEff(SES_pd_mod, R = R_sample, ran.eff = "Realm")
eff_pd_ses <- semEff(boot_mod_pd_ses)

Sys.time()


################## Save outputs for figures ################


# fd
saveRDS(fric_mod, "Output/10_SEM_FD_realm_model.rds")
saveRDS(eff_fd, "Output/10_SEM_FD_realm_boot_effects.rds")
# pd
saveRDS(pd_mod, "Output/10_SEM_PD_realm_model.rds")
saveRDS(eff_pd, "Output/10_SEM_PD_realm_boot_effects.rds")

# SES_fd
saveRDS(SES_fric_mod, "Output/10_SEM_sesFD_realm_model.rds")
saveRDS(eff_fd_ses, "Output/10_SEM_sesFD_realm_boot_effects.rds")
# SES_pd
saveRDS(SES_pd_mod, "Output/10_SEM_sesPD_realm_model.rds")
saveRDS(eff_pd_ses, "Output/10_SEM_sesPD_realm_boot_effects.rds")


