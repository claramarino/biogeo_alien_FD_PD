# Compute SEM Models for 407 world oceanic islands 
# with at least 4 alien bird species
# remove human variables


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
colSums(is.na(all_ok))

# remove CP and TFSI for main models because of NA
all_ok <- all_ok %>% select(-c(ARCHIP, ISLAND, CP, TSFI))

###################### Model specification #########################

# Specify submodels with realm as random effect
colnames(all_ok)


# No human variable

# biotic context depends on biogeographic variables
# native SR drives native fd & PD
nat_SR <- lmer(SR_nat ~ Area + Dist + Lat + Elev + SLMP + 
                 (1|Realm), data = all_ok)
nat_PD <- lmer(PD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                 (1|Realm), data = all_ok)
nat_FD <- lmer(FD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat +
                   (1|Realm), data = all_ok)
# alien diversities depend on all contextual variables
alien_SR <- lmer(SR_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                   (1|Realm), data = all_ok)
alien_PD <- lmer(PD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + PD_nat +
                   SR_alien + 
                   (1|Realm), data = all_ok)
alien_FD <- lmer(FD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + FD_nat +
                     SR_alien +
                     (1|Realm), data = all_ok)

# build total model
# one for FD and one for pd

pd_mod <- psem(nat_SR, nat_PD, alien_SR, 
               alien_PD, data = all_ok)

basisSet(pd_mod)
d1 <- dSep(pd_mod)
fisherC(pd_mod)
d1[which.max(abs(d1$Crit.Value)),] 
# adding correlated error between native diversities and human factors
# until the pval is above .05 => ok here

fd_mod <- psem(nat_SR, nat_FD, alien_SR, 
                 alien_FD, data = all_ok)
d1 <- dSep(fd_mod)
fisherC(fd_mod) 
d1[which.max(abs(d1$Crit.Value)),] 
# adding correlated error between native diversities and human factors
# until the pval is above .05 => ok here


############## Simple model outputs #############

# for PD 
rsquared(pd_mod) # Fort effet sur la richesse spé native et alien cf. R² marginal/conditional
coefs <- coefs(pd_mod)
coefs[coefs$P.Value<0.05,]


# for FD 
rsquared(fd_mod) # Mm chose qu'avant + effet assez important sur FD (R²+0.2)
coefs <- coefs(fd_mod) 
coefs[coefs$P.Value<0.05,]


########### Response variables: SES-PD and SES-FD ###################

# Hypotheses are the same as above but we take the SES as response variables
# and not the diversity values per se

# no human variables

# biotic context depends on biogeographic variables
# native SR drives native FD & PD
nat_SR <- lmer(SR_nat ~ Area + Dist + Lat + Elev + SLMP + 
                 (1|Realm), data = all_ok)
nat_SES_PD <- lmer(SES_PD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                     (1|Realm), data = all_ok)
nat_SES_FD <- lmer(SES_FD_nat ~ Area + Dist + Lat + Elev + SLMP + SR_nat +
                       (1|Realm), data = all_ok)
# alien diversities depend on all contextual variables
alien_SR <- lmer(SR_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat +
                   (1|Realm), data = all_ok)
alien_SES_PD <- lmer(SES_PD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                       SES_PD_nat + SR_alien + 
                       (1|Realm), data = all_ok)
alien_SES_FD <- lmer(SES_FD_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + 
                       SES_FD_nat + SR_alien +
                         (1|Realm), data = all_ok)

# build total model
# one for SES_PD and one for SES_FD
# alien diversities depend on all contextual variables
alien_SR2 <- lmer(SR_alien ~ Area + Dist + Lat + Elev + SLMP + SR_nat + SES_PD_nat +
                    (1|Realm), data = all_ok)
SES_pd_mod <- psem(nat_SR, nat_SES_PD, alien_SR, alien_SES_PD,
                   #SES_PD_nat %~~% SR_alien,
                   data = all_ok)

d1 <- dSep(SES_pd_mod)
fisherC(SES_pd_mod) # F.C = 40.947, df = 34, P = 0.192 
d1[which.max(abs(d1$Crit.Value)),] 
# adding correlated error between native diversities and human factors
# until the pval is above .05
# + SES PD nat as a predictor of alien SR

SES_fd_mod <- psem(nat_SR, nat_SES_FD, alien_SR, alien_SES_FD,
                     data = all_ok)
d1 <- dSep(SES_fd_mod)
fisherC(SES_fd_mod) # F.C = 50.754, df = 40, P = 0.119
d1[which.max(abs(d1$Crit.Value)),] 
# adding correlated error between native diversities and human factors
# until the pval is above .05 (less additions than for FD)



############## Simple model outputs #############

# for PD 
rsquared(SES_pd_mod) # Fort effet sur la richesse spé native et alien cf. R² marginal/conditional
coefs <- coefs(SES_pd_mod)
coefs[coefs$P.Value<0.05,]


# for FD 
rsquared(SES_fd_mod) # Mm chose qu'avant + effet assez important sur FD (R²+0.2)
coefs <- coefs(SES_fd_mod) 
coefs[coefs$P.Value<0.05,]


########## Get direct/ indirect effects of the 4 models #############

#install.packages("semEff")

# resample for confidence intervals
R_sample = 10000

Sys.time()

boot_mod_fd <- bootEff(fd_mod, R = R_sample, ran.eff = "Realm")
eff_fd <- semEff(boot_mod_fd)

Sys.time()

boot_mod_pd <- bootEff(pd_mod, R = R_sample, ran.eff = "Realm")
eff_pd <- semEff(boot_mod_pd)

Sys.time()

boot_mod_fd_ses <- bootEff(SES_fd_mod, R = R_sample, ran.eff = "Realm")
eff_fd_ses <- semEff(boot_mod_fd_ses)

Sys.time()

boot_mod_pd_ses <- bootEff(SES_pd_mod, R = R_sample, ran.eff = "Realm")
eff_pd_ses <- semEff(boot_mod_pd_ses)

Sys.time()


################## Save outputs for figures ################


# fd
saveRDS(fd_mod, "Output/13_SEM_FD_realm_model_no_human.rds")
saveRDS(eff_fd, "Output/13_SEM_FD_realm_boot_effects_no_human.rds")
# pd
saveRDS(pd_mod, "Output/13_SEM_PD_realm_model_no_human.rds")
saveRDS(eff_pd, "Output/13_SEM_PD_realm_boot_effects_no_human.rds")

# SES_fd
saveRDS(SES_fd_mod, "Output/13_SEM_sesFD_realm_model_no_human.rds")
saveRDS(eff_fd_ses, "Output/13_SEM_sesFD_realm_boot_effects_no_human.rds")
# SES_pd
saveRDS(SES_pd_mod, "Output/13_SEM_sesPD_realm_model_no_human.rds")
saveRDS(eff_pd_ses, "Output/13_SEM_sesPD_realm_boot_effects_no_human.rds")


#################### Plot Figure S4 #######################


# load models
# fd
mod_fd <- readRDS("Output/13_SEM_FD_realm_model_no_human.rds")
eff_fd <- readRDS("Output/13_SEM_FD_realm_boot_effects_no_human.rds")
# pd
mod_pd <- readRDS("Output/13_SEM_PD_realm_model_no_human.rds")
eff_pd <- readRDS("Output/13_SEM_PD_realm_boot_effects_no_human.rds")

# SES_fd
mod_fd_ses <- readRDS("Output/13_SEM_sesFD_realm_model_no_human.rds")
eff_fd_ses <- readRDS("Output/13_SEM_sesFD_realm_boot_effects_no_human.rds")
# SES_pd
mod_pd_ses <- readRDS("Output/13_SEM_sesPD_realm_model_no_human.rds")
eff_pd_ses <- readRDS("Output/13_SEM_sesPD_realm_boot_effects_no_human.rds")


# Extract effects as data.frame for figures

t_eff_fd <- as.data.frame(eff_fd$Summary$FD.alien)
t_eff_pd <- as.data.frame(eff_pd$Summary$PD.alien)
t_eff_fd_ses <- as.data.frame(eff_fd_ses$Summary$SES.FD.alien)
t_eff_pd_ses <- as.data.frame(eff_pd_ses$Summary$SES.PD.alien)



# use a function to reshape tables
shape_tbl <- function(t_pourri){
  t_mieux <- t_pourri[-c(1, 10,17,26),-c(3,5,7,9,12)]
  
  names(t_mieux) <- c("Type", "Var", "Effect", "Bias", "Std_Err", "Lower_CI", "Upper_CI", "signif")
  t_encoremieux <- t_mieux %>% 
    mutate(Type = c(rep("Direct", 8), rep("Indirect", 6),
                    rep("Total",8), rep("Mediator", 3))) %>%
    mutate_at(c("Effect", "Bias", "Std_Err", "Lower_CI", "Upper_CI"),
              as.numeric) %>%
    mutate_at(c("Type", "Var", "signif"), as.character) %>%
    mutate(Var = gsub(" ", "", Var))%>%
    mutate(Var = factor(Var, levels = c(
      "Lat", "Area", "Dist", "SLMP", "Elev", "GDP", "pop", "connect" ,"modif.change", 
      "static.modif", "SR.nat", "PD.nat","SES.PD.nat","FD.nat", "SES.FD.nat", 
      "SR.alien"), ordered = T))%>%
    mutate(signif = factor(signif, levels = c("*", ""), ordered = T))
  return(t_encoremieux)
}


fd_fig <- shape_tbl(t_eff_fd) %>% mutate(response = "FD")
pd_fig <- shape_tbl(t_eff_pd) %>% mutate(response = "PD")
ses_fd_fig <- shape_tbl(t_eff_fd_ses) %>% mutate(response = "SES_FD")
ses_pd_fig <- shape_tbl(t_eff_pd_ses) %>% mutate(response = "SES_PD")


# FD and SES-FD
fd_and_ses <- bind_rows(fd_fig, ses_fd_fig) %>%
  filter(Type %in% c("Direct", "Indirect", "Total"))

# PD and SES-PD
pd_and_ses <- bind_rows(pd_fig, ses_pd_fig) %>%
  filter(Type %in% c("Direct", "Indirect", "Total"))

pd = position_dodge(.5)

p1bis <- ggplot(data = fd_and_ses, aes(x = Var, y = Effect, color = response,
                                       alpha = signif))+
  geom_abline(slope = 0, intercept = 0, linetype = 2)+
  geom_errorbar(aes(ymin=Lower_CI, ymax=Upper_CI),
                width=.4, position = pd, size = .6) +
  geom_point(position = pd, size=2) +
  scale_alpha_manual(values = c(1, .3))+
  scale_color_manual(values=c('turquoise3', 'indianred3'))+
  facet_wrap(~Type)+
  coord_flip()+
  theme_light()+
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12))
p2bis <- ggplot(data = pd_and_ses, aes(x = Var, y = Effect, color = response, 
                                       alpha = signif))+
  geom_abline(slope = 0, intercept = 0, linetype = 2)+
  geom_errorbar(aes(ymin=Lower_CI, ymax=Upper_CI),
                width=.4, position = pd, size = .6) +
  geom_point(position = pd, size=2) +
  scale_alpha_manual(values = c(1, .3))+
  scale_color_manual(values=c('springgreen3', 'violetred3'))+
  facet_wrap(~Type)+
  coord_flip()+
  theme_light()+
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12))

pdf("Fig/13_Effects_dir_indir_tot_fd_no_human.pdf", 6, 3)
p1bis
dev.off()
pdf("Fig/13_Effects_dir_indir_tot_pd_no_human.pdf", 6, 3)
p2bis
dev.off()

