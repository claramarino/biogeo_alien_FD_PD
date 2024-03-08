# SEM figure & tables
rm(list=ls())

library(semEff)
library(tidyverse)


# load models

# save objects for figures 
# fd
mod_fd <- readRDS("Output/01_SEM_Fric_realm_model")
eff_fd <- readRDS("Output/01_SEM_Fric_realm_boot_effects")
# pd
mod_pd <- readRDS("Output/01_SEM_PD_realm_model")
eff_pd <- readRDS("Output/01_SEM_PD_realm_boot_effects")

# SES_fd
mod_fd_ses <- readRDS("Output/01_SEM_sesFric_realm_model")
eff_fd_ses <- readRDS("Output/01_SEM_sesFric_realm_boot_effects")
# SES_pd
mod_pd_ses <- readRDS("Output/01_SEM_sesPD_realm_model")
eff_pd_ses <- readRDS("Output/01_SEM_sesPD_realm_boot_effects")


eff_fd$Summary

# Extract effects as data.frame for figures

t_eff_fd <- as.data.frame(eff_fd$Summary$fric.alien)
t_eff_pd <- as.data.frame(eff_pd$Summary$PD.alien)
t_eff_fd_ses <- as.data.frame(eff_fd_ses$Summary$SES.fric.alien)
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
      "static.modif", "SR.nat", "PD.nat","SES.PD.nat","fric.nat", "SES.fric.nat", 
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

# fric and ses fric


fd_and_ses <- bind_rows(fd_fig, ses_fd_fig) %>%
  filter(Type %in% c("Direct", "Indirect", "Total"))
                      

pd = position_dodge(.5)

str(fd_and_ses)

p1 <- ggplot(data = fd_and_ses, aes(x = Var, y = Effect, color = response, 
                              linetype = signif, shape = signif))+
  geom_abline(slope = 0, intercept = 0, linetype = 2)+
  geom_errorbar(aes(ymin=Lower_CI, ymax=Upper_CI),
                width=.4, position = pd, linewidth = .6) +
  geom_point(position = pd, size=2) +
  scale_shape_manual(values=c(16, 1))+
  scale_color_manual(values=c('turquoise3', 'indianred3'))+
  facet_wrap(~Type)+
  coord_flip()+
  theme_light()+
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12))

p1


pd_and_ses <- bind_rows(pd_fig, ses_pd_fig) %>%
  filter(Type %in% c("Direct", "Indirect", "Total"))

p2 <- ggplot(data = pd_and_ses, aes(x = Var, y = Effect, color = response, 
                                    linetype = signif, shape = signif))+
  geom_abline(slope = 0, intercept = 0, linetype = 2)+
  geom_errorbar(aes(ymin=Lower_CI, ymax=Upper_CI),
                width=.4, position = pd, size = .6) +
  geom_point(position = pd, size=2) +
  scale_shape_manual(values=c(16, 1))+
  scale_color_manual(values=c('springgreen3', 'violetred3'))+
  facet_wrap(~Type)+
  coord_flip()+
  theme_light()+
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12))


p2

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

p2bis



pdf("Fig/16_Effects_dir_indir_tot_fd.pdf", 8, 4)
p1bis
dev.off()
pdf("Fig/16_Effects_dir_indir_tot_pd.pdf", 8, 4)
p2bis
dev.off()




p1ter <- ggplot(data = fd_and_ses, aes(x = Var, y = Effect, color = response))+
  geom_abline(slope = 0, intercept = 0, linetype = 2)+
  geom_errorbar(aes(ymin=Lower_CI, ymax=Upper_CI),
                width=.4, position = pd, size = .6) +
  geom_point(position = pd, size=2) +
  scale_color_manual(values=c('turquoise3', 'indianred3'))+
  facet_wrap(~Type)+
  coord_flip()+
  theme_light()+
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12))
p2ter <- ggplot(data = pd_and_ses, aes(x = Var, y = Effect, color = response))+
  geom_abline(slope = 0, intercept = 0, linetype = 2)+
  geom_errorbar(aes(ymin=Lower_CI, ymax=Upper_CI),
                width=.4, position = pd, size = .6) +
  geom_point(position = pd, size=2) +
  scale_color_manual(values=c('springgreen3', 'violetred3'))+
  facet_wrap(~Type)+
  coord_flip()+
  theme_light()+
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12))


pdf("Fig/31_Effects_dir_indir_tot_fd.pdf", 8, 4)
p1ter
dev.off()
pdf("Fig/31_Effects_dir_indir_tot_pd.pdf", 8, 4)
p2ter
dev.off()




fisherC(random_mod_fd)

rsquared(random_mod_fd)$Response
round(rsquared(random_mod_fd)$Marginal,2)
round(rsquared(random_mod_fd)$Conditional,2)

coefs <- coefs(random_mod_fd) 
coefs[coefs$P.Value<0.05,]


Dir_effect <- getDirEff(eff_fd, "fric") # Extract direct effects.
fric_dir <- as.data.frame(Dir_effect) %>%
  rownames_to_column("Var") %>%
  rename(Effect = Dir_effect) %>%
  mutate(type="Direct", 
         response = "FD")

Indir_effect <- getIndEff(eff_fd, "fric") # Extract indirect effects.
fric_indir <- as.data.frame(Indir_effect) %>%
  rownames_to_column("Var") %>%
  rename(Effect = Indir_effect) %>%
  mutate(type="Indirect", 
         response = "FD")

Effect = getTotEff(eff_fd, "fric") # Extract total effects.
fric_tot <- as.data.frame(Effect) %>%
  rownames_to_column("Var") %>%
  mutate(type="Total", 
         response = "FD")
# calculate relative direct effect
rel = sum(abs(fric_tot$Effect))
fric_tot_rel <- fric_tot %>%
  mutate(Effect = abs(Effect)/rel) %>%
  mutate(type="TotalRel")


fric_all <- bind_rows(fric_dir, fric_indir, fric_tot)

fric_all$Var <- factor(fric_all$Var, levels = fric_all$Var[order(fric_all$Effect, decreasing = TRUE)])

ggplot(data = fric_all)+
  geom_bar(aes(x = Var, y = Effect), 
           stat = "identity", position = "dodge")+
  facet_wrap(~type)+
  coord_flip()


library(piecewiseSEM)
rsquared(mod_fd)
rsquared(mod_pd)
rsquared(random_mod_pd)
rsquared(random_mod_pd_ses) # Mm chose qu'avant + effet assez important sur Fric (R²+0.1)

coefs <- coefs(random_mod_fd) 
coefs[coefs$P.Value<0.05,]

eff_fd$Summary$fric


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


eff_fd$Summary$fric
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
