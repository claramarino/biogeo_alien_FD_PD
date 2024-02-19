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

###################### Model specification #########################

colnames(all_ok)
colSums(is.na(all_ok))


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
                   nat_SR, nat_PD, alien_SR, alien_PD,
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
                 nat_SR, nat_fric, alien_SR, alien_fric,
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
saveRDS(fric_mod, "Output/15_SEM_Fric_realm_model")
saveRDS(eff_fd, "Output/15_SEM_Fric_realm_boot_effects")
# pd
saveRDS(pd_mod, "Output/15_SEM_PD_realm_model")
saveRDS(eff_pd, "Output/15_SEM_PD_realm_boot_effects")

# SES_fd
saveRDS(SES_fric_mod, "Output/15_SEM_sesFric_realm_model")
saveRDS(eff_fd_ses, "Output/15_SEM_sesFric_realm_boot_effects")
# SES_pd
saveRDS(SES_pd_mod, "Output/15_SEM_sesPD_realm_model")
saveRDS(eff_pd_ses, "Output/15_SEM_sesPD_realm_boot_effects")


################## Compare diversities alien / native ################

summary(all)

all_lg <- all %>% 
  rename(SESPD_alien = SES_PD_alien,
         SESFD_alien = SES_FD_alien,
         SESPD_nat = SES_PD_nat,
         SESFD_nat = SES_FD_nat) %>%
  pivot_longer(cols = c(SR_nat, SR_alien:SESFD_nat), 
               names_to = "metric", values_to = "value") %>%
  separate(metric, c("metric","community"))


ggplot(data = all_lg)+
  geom_boxplot(aes(x = community, y = value))+
  facet_wrap(~metric, scales = "free_y")

ggplot(all, aes(x=PD_alien, y = PD_nat))+
  geom_point()+ geom_smooth(method = "lm")


ggplot(all, aes(x=FD_alien, y = FD_nat))+
  geom_point()+ geom_smooth(method = "lm")


################ Add colonization pressure variable ################

# For a subset of islands (n=96), we could obtain data regarding
# the colonization pressure (CP), that is the total number of species
# introduced to an area (Blackburn et al., 2020)

# We thus perform SEM on those 96 islands with and without the CP variable 
# to check how this factor can influence our results.

# we still keep the total SEMs on the 407 islands without the CP XXXXXX


# Filter the islands with CP data
isl_CP <- all_ok %>% filter(!is.na(CP))
table(isl_CP$Realm)
table(all_ok$Realm)


# human variables depend on biogeographic var
# and are related to each other
pop <- lmer(pop ~ Area + Lat + SLMP + Dist +
              (1|Realm), data = isl_CP)
GDP <- lmer(GDP ~ Area + pop + Dist + Elev +
              (1|Realm), data = isl_CP)
connect <- lmer(connect ~ Area + Lat + SLMP + Dist + pop + 
                  (1|Realm), data = isl_CP)
static <- lmer(static_modif ~ GDP + pop + modif_change + Area + Lat + connect + 
                 (1|Realm), data = isl_CP)
modif <- lmer(modif_change ~ GDP + pop + Elev  + SLMP+ connect +
                (1|Realm), data = isl_CP)
# biotic context depends on biogeographic variables
# native SR drives native fric & PD
nat_SR <- lmer(SR_nat ~ Area + Dist + Lat + Elev + SLMP + 
                 (1|Realm), data = isl_CP)
nat_PD <- lmer(PD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                 (1|Realm), data = isl_CP)
nat_fric <- lmer(FD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat +
                   (1|Realm), data = isl_CP)
# alien diversities depend on all contextual variables
alien_SR <- lmer(SR_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                   pop + GDP + modif_change + static_modif + connect +
                   (1|Realm), data = isl_CP)
alien_PD <- lmer(PD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + PD_nat +
                   pop + GDP + modif_change + static_modif + connect + SR_alien + 
                   (1|Realm), data = isl_CP)
alien_fric <- lmer(FD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + FD_nat +
                     pop + GDP + modif_change + static_modif + connect + SR_alien +
                     (1|Realm), data = isl_CP)

# build total model
# one for fric and one for pd

pd_mod <- psem(pop, GDP, connect, static, modif, 
               nat_SR, nat_PD, alien_SR, alien_PD,
               SR_nat %~~% connect,
               PD_nat %~~% pop,
               PD_nat %~~% connect,
               PD_nat %~~% GDP,
               data = isl_CP)

basisSet(pd_mod)
d1 <- dSep(pd_mod)
fisherC(pd_mod) # F.C = 44.366, df = 36, P = 0.16 
d1[which.max(abs(d1$Crit.Value)),] 
# adding correlated error between native diversities and human factors
# until the pval is above .05

fric_mod <- psem(pop, GDP, connect, static, modif, 
                 nat_SR, nat_fric, alien_SR, alien_fric,
                 SR_nat %~~% connect,
                 FD_nat %~~% pop,
                 FD_nat %~~% static_modif, 
                 FD_nat %~~% connect,
                 data = isl_CP)
d1 <- dSep(fric_mod)
fisherC(fric_mod) # F.C = 42.253, df = 36, P = 0.219 
d1[which.max(abs(d1$Crit.Value)),] 
# adding correlated error between native diversities and human factors
# until the pval is above .05


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
              (1|Realm), data = isl_CP)
GDP <- lmer(GDP ~ Area + pop + Dist + Elev +
              (1|Realm), data = isl_CP)
connect <- lmer(connect ~ Area + Lat + SLMP + Dist + pop + 
                  (1|Realm), data = isl_CP)
static <- lmer(static_modif ~ GDP + pop + modif_change + Area + Lat + connect + 
                 (1|Realm), data = isl_CP)
modif <- lmer(modif_change ~ GDP + pop + Elev  + SLMP+ connect +
                (1|Realm), data = isl_CP)
# biotic context depends on biogeographic variables
# native SR drives native fric & PD
nat_SR <- lmer(SR_nat ~ Area + Dist + Lat + Elev + SLMP + 
                 (1|Realm), data = isl_CP)
nat_SES_PD <- lmer(SES_PD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                     (1|Realm), data = isl_CP)
nat_SES_fric <- lmer(SES_FD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat +
                       (1|Realm), data = isl_CP)
# alien diversities depend on all contextual variables
alien_SR <- lmer(SR_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                   pop + GDP + modif_change + static_modif + connect +
                   (1|Realm), data = isl_CP)
alien_SES_PD <- lmer(SES_PD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + SES_PD_nat +
                       pop + GDP + modif_change + static_modif + connect + SR_alien + 
                       (1|Realm), data = isl_CP)
alien_SES_fric <- lmer(SES_FD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + SES_FD_nat +
                         pop + GDP + modif_change + static_modif + connect + SR_alien +
                         (1|Realm), data = isl_CP)

# build total model
# one for SES_PD and one for SES_fric
# alien diversities depend on all contextual variables
alien_SR2 <- lmer(SR_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + SES_PD_nat +
                    pop + GDP + modif_change + static_modif + connect +
                    (1|Realm), data = isl_CP)
SES_pd_mod <- psem(pop, GDP, connect, static, modif, 
                   nat_SR, nat_SES_PD, alien_SR2, alien_SES_PD,
                   SR_nat %~~% connect,
                   SES_PD_nat %~~% pop,
                   SES_PD_nat %~~% connect,
                   SES_PD_nat %~~% GDP,
                   data = isl_CP)

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
                     data = isl_CP)
d1 <- dSep(SES_fric_mod)
fisherC(SES_fric_mod) # F.C = 51.509, df = 40, P = 0.105 
d1[which.max(abs(d1$Crit.Value)),] 
# adding correlated error between native diversities and human factors
# until the pval is above .05 (less additions than for Fric)


# carreful: pbm with singularity of models
# check the effects of some var?
isSingular(alien_SES_PD)
isSingular(alien_SES_fric)

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
saveRDS(fric_mod, "Output/01_SEM_R1_CP_FD_model")
saveRDS(eff_fd, "Output/01_SEM_R1_CP_FD_boot_eff")
# pd
saveRDS(pd_mod, "Output/01_SEM_R1_CP_PD_model")
saveRDS(eff_pd, "Output/01_SEM_R1_CP_PD_boot_eff")

# SES_fd
saveRDS(SES_fric_mod, "Output/01_SEM_R1_CP_SESFD_model")
saveRDS(eff_fd_ses, "Output/01_SEM_R1_CP_SESFD_boot_eff")
# SES_pd
saveRDS(SES_pd_mod, "Output/01_SEM_R1_CP_SESPD_model")
saveRDS(eff_pd_ses, "Output/01_SEM_R1_CP_SESPD_boot_eff")




