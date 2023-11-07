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

# contextual variables
var <- readRDS("Data/07_oce_isl.rds") 

# alien diversities
fric_alien <- readRDS("Data/13_func_div_ses_world.rds")
pd_alien <- readRDS("Data/11_pd_by_isl_world.rds")

# native diversities
fric_nat <- readRDS("Data/13_fric_+_ses_natives_isl_over3.rds")
pd_nat <- readRDS("Data/11_pd_+_ses_natives_world.rds")


######### Clean and join tables ############


var_ok <- var %>% 
  mutate(connect = nb_ap + nb_ap_buff_100) %>% # transform connect as slmp 
  # remove var not used in models
  select(-c(varP, varT, HFI, GDP_mean, GDP_max, GDP_0_free, intact, Long,
            GMMC, nb_ap, nb_ap_buff_50, nb_ap_buff_100, dist_ap)) %>%
  rename(SR_nat = native_sp_rich,
         GDP = GDP_sum) %>% 
  mutate(ID = as.character(ID))


nat_all <- left_join(pd_nat, fric_nat) %>%
  rename(PD_nat = PD_mean, 
         SES_PD_nat = ses_PD, 
         ID = isl) %>%
  select(-SR)

alien_all <- left_join(fric_alien, pd_alien %>% select(-c(SR))) %>%
  select(isl:fric, SES_fric, PD_mean, ses_PD) %>%
  rename(fric_alien = fric, 
         SES_fric_alien = SES_fric, 
         PD_alien = PD_mean,
         SES_PD_alien = ses_PD, 
         SR_alien = SR,
         ID = isl)

all = left_join(var_ok, left_join(alien_all, nat_all))
  
# scale variables after log transformation
colnames(all)
all_ok <- all %>%
  mutate_at(c("Area", "Dist", "fric_alien", "fric_nat", "PD_alien", 
              "PD_nat", "SR_alien", "SR_nat"), log) %>%
  mutate_at(c("pop", "Elev","GDP", "connect"), ~log(.+1)) %>%
  mutate(Lat = abs(Lat)) %>%
  # scale variables before SEM
  mutate_if(is.numeric, ~as.numeric(scale(.)))


table(all_ok$Realm)

colSums(is.na(all_ok))

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
nat_fric <- lmer(fric_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat +
                   (1|Realm), data = all_ok)
# alien diversities depend on all contextual variables
alien_SR <- lmer(SR_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                   pop + GDP + modif_change + static_modif + connect +
                   (1|Realm), data = all_ok)
alien_PD <- lmer(PD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + PD_nat +
                   pop + GDP + modif_change + static_modif + connect + SR_alien + 
                   (1|Realm), data = all_ok)
alien_fric <- lmer(fric_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + fric_nat +
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
fisherC(pd_mod) # F.C = 44.366, df = 36, P = 0.16 
d1[which.max(abs(d1$Crit.Value)),] 
# adding correlated error between native diversities and human factors
# until the pval is above .05

fric_mod <- psem(pop, GDP, connect, static, modif, 
                 nat_SR, nat_fric, alien_SR, alien_fric,
                 SR_nat %~~% connect,
                 fric_nat %~~% pop,
                 fric_nat %~~% static_modif, 
                 fric_nat %~~% connect,
                 data = all_ok)
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


########### Random effect on SES pd/fric ###################

# Starting with final total model + adding random effect of Realms everywhere
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
nat_SES_fric <- lmer(SES_fric_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat +
                   (1|Realm), data = all_ok)
# alien diversities depend on all contextual variables
alien_SR <- lmer(SR_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                   pop + GDP + modif_change + static_modif + connect +
                   (1|Realm), data = all_ok)
alien_SES_PD <- lmer(SES_PD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + SES_PD_nat +
                   pop + GDP + modif_change + static_modif + connect + SR_alien + 
                   (1|Realm), data = all_ok)
alien_SES_fric <- lmer(SES_fric_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + SES_fric_nat +
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
                 SES_fric_nat %~~% pop,
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
         SESfric_alien = SES_fric_alien,
         SESPD_nat = SES_PD_nat,
         SESfric_nat = SES_fric_nat) %>%
  pivot_longer(cols = c(SR_nat, SR_alien:SESfric_nat), 
               names_to = "metric", values_to = "value") %>%
  separate(metric, c("metric","community"))


ggplot(data = all_lg)+
  geom_boxplot(aes(x = community, y = value))+
  facet_wrap(~metric, scales = "free_y")

ggplot(all, aes(x=PD_alien, y = PD_nat))+
  geom_point()+ geom_smooth(method = "lm")


ggplot(all, aes(x=fric_alien, y = fric_nat))+
  geom_point()+ geom_smooth(method = "lm")
