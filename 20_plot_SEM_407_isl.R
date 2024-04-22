# SEM figure & tables
rm(list=ls())

library(semEff)
library(tidyverse)


# load models
# fd
mod_fd <- readRDS("Output/10_SEM_FD_realm_model.rds")
eff_fd <- readRDS("Output/10_SEM_FD_realm_boot_effects.rds")
# pd
mod_pd <- readRDS("Output/10_SEM_PD_realm_model.rds")
eff_pd <- readRDS("Output/10_SEM_PD_realm_boot_effects.rds")

# SES_fd
mod_fd_ses <- readRDS("Output/10_SEM_sesFD_realm_model.rds")
eff_fd_ses <- readRDS("Output/10_SEM_sesFD_realm_boot_effects.rds")
# SES_pd
mod_pd_ses <- readRDS("Output/10_SEM_sesPD_realm_model.rds")
eff_pd_ses <- readRDS("Output/10_SEM_sesPD_realm_boot_effects.rds")

eff_fd$Summary

# Extract effects as data.frame for figures

t_eff_fd <- as.data.frame(eff_fd$Summary$FD.alien)
t_eff_pd <- as.data.frame(eff_pd$Summary$PD.alien)
t_eff_fd_ses <- as.data.frame(eff_fd_ses$Summary$SES.FD.alien)
t_eff_pd_ses <- as.data.frame(eff_pd_ses$Summary$SES.PD.alien)

# use a function to reshape tables
shape_tbl <- function(t_pourri){
  if (nrow(t_pourri)==49) {i=0} else {i=1}
  t_mieux <- t_pourri[-c(1, 15,27+i,41+i),-c(3,5,7,9,12)]
  
  names(t_mieux) <- c("Type", "Var", "Effect", "Bias", "Std_Err", "Lower_CI", "Upper_CI", "signif")
  t_encoremieux <- t_mieux %>% 
    mutate(Type = c(rep("Direct", 13), rep("Indirect", 11+i),
                    rep("Total",13), rep("Mediator", 8))) %>%
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

eff_pd$Summary$SR.alien
eff_pd$Summary$PD.nat

# combine tables for ggplot

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

pdf("Fig/20_Effects_dir_indir_tot_fd.pdf", 8, 4)
p1bis
dev.off()
pdf("Fig/20_Effects_dir_indir_tot_pd.pdf", 8, 4)
p2bis
dev.off()




#######################################################################

Dir_effect <- getDirEff(eff_fd, "FD") # Extract direct effects.
FD_dir <- as.data.frame(Dir_effect) %>%
  rownames_to_column("Var") %>%
  rename(Effect = Dir_effect) %>%
  mutate(type="Direct", 
         response = "FD")

Indir_effect <- getIndEff(eff_fd, "FD") # Extract indirect effects.
FD_indir <- as.data.frame(Indir_effect) %>%
  rownames_to_column("Var") %>%
  rename(Effect = Indir_effect) %>%
  mutate(type="Indirect", 
         response = "FD")

Effect = getTotEff(eff_fd, "FD") # Extract total effects.
FD_tot <- as.data.frame(Effect) %>%
  rownames_to_column("Var") %>%
  mutate(type="Total", 
         response = "FD")
# calculate relative direct effect
rel = sum(abs(FD_tot$Effect))
FD_tot_rel <- FD_tot %>%
  mutate(Effect = abs(Effect)/rel) %>%
  mutate(type="TotalRel")


FD_all <- bind_rows(FD_dir, FD_indir, FD_tot)

FD_all$Var <- factor(FD_all$Var, levels = FD_all$Var[order(FD_all$Effect, decreasing = TRUE)])

ggplot(data = FD_all)+
  geom_bar(aes(x = Var, y = Effect), 
           stat = "identity", position = "dodge")+
  facet_wrap(~type)+
  coord_flip()



# fig for presentation
# remove SES
pdf("Fig/31_Pres_SFE2_Effects_fd.pdf", 8, 4)
ggplot(data = fd_fig %>%
         filter(Type %in% c("Direct", "Indirect", "Total")), 
       aes(x = Var, y = Effect, alpha = signif))+
  geom_abline(slope = 0, intercept = 0, linetype = 2)+
  geom_errorbar(aes(ymin=Lower_CI, ymax=Upper_CI),
                width=.4, size = .6, color = "violetred3") +
  geom_point(size=2, color = "violetred3") +
  scale_alpha_manual(values = c(1, .3))+
  facet_wrap(~Type)+
  coord_flip()+
  theme_light()+
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12))
dev.off()




# draw the path
library(piecewiseSEM)
rsquared(random_mod_fd)

rsquared(random_mod_pd)

coef <- coefs(random_mod_fd)

coef[coef$P.Value<0.05,]


eff_fd$Summary$FD
eff_pd$Summary$PD.mean

eff_fd$Summary$SR
eff_fd$Summary$native.sp.rich
eff_fd$Summary$modif.change
eff_fd$Summary$static.modif
eff_fd$Summary$connect
eff_fd$Summary$GDP
eff_fd$Summary$pop

eff_pd$Summary$SR
eff_pd$Summary$native.sp.rich
eff_pd$Summary$modif.change
eff_pd$Summary$static.modif
eff_pd$Summary$connect
eff_pd$Summary$GDP
eff_pd$Summary$pop
