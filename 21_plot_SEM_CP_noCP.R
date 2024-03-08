# SEM figure & tables
# with and without CP variables
# on the subset of 96 islands

rm(list=ls())

library(semEff)
library(tidyverse)


# load models

# FD without CP
mod_fd <- readRDS("Output/11_SEM_noCP_FD_model.rds")
eff_fd <- readRDS("Output/11_SEM_noCP_FD_boot_eff.rds")
# PD without CP
mod_pd <- readRDS("Output/11_SEM_noCP_PD_model.rds")
eff_pd <- readRDS("Output/11_SEM_noCP_PD_boot_eff.rds")

# FD with CP
mod_fd_CP <- readRDS("Output/11_SEM_CP_FD_model.rds")
eff_fd_CP <- readRDS("Output/11_SEM_CP_FD_boot_eff.rds")
# PD with CP
mod_pd_CP <- readRDS("Output/11_SEM_CP_PD_model.rds")
eff_pd_CP <- readRDS("Output/11_SEM_CP_PD_boot_eff.rds")


eff_fd$Summary

# Extract effects as data.frame for figures

t_eff_fd <- as.data.frame(eff_fd$Summary$FD.alien)
t_eff_pd <- as.data.frame(eff_pd$Summary$PD.alien)
t_eff_fd_cp <- as.data.frame(eff_fd_CP$Summary$FD.alien)
t_eff_pd_cp <- as.data.frame(eff_pd_CP$Summary$PD.alien)

# use a function to reshape tables
shape_tbl <- function(t_pourri){
  t_mieux <- t_pourri[-c(1,15,27,41),-c(3,5,7,9,12)]
  
  names(t_mieux) <- c("Type", "Var", "Effect", "Bias", "Std_Err", "Lower_CI", "Upper_CI", "signif")
  t_encoremieux <- t_mieux %>% 
    mutate(Type = c(rep("Direct", 13), rep("Indirect", 11),
                    rep("Total",13), rep("Mediator", 8))) %>%
    mutate_at(c("Effect", "Bias", "Std_Err", "Lower_CI", "Upper_CI"),
              as.numeric) %>%
    mutate_at(c("Type", "Var", "signif"), as.character) %>%
    mutate(Var = gsub(" ", "", Var))%>%
    mutate(Var = factor(Var, levels = c(
      "Lat", "Area", "Dist", "SLMP", "Elev", "GDP", "pop", "connect" ,"modif.change", 
      "static.modif", "CP","SR.nat", "PD.nat","SES.PD.nat","FD.nat", "SES.FD.nat", 
      "SR.alien"), ordered = T))%>%
    mutate(signif = factor(signif, levels = c("*", ""), ordered = T))
  return(t_encoremieux)

}

# use a function to reshape tables with CP
shape_tbl_CP <- function(t_pourri){
  t_mieux <- t_pourri[-c(1,16,29,44), -c(3,5,7,9,12)]
  
  names(t_mieux) <- c("Type", "Var", "Effect", "Bias", "Std_Err", "Lower_CI", "Upper_CI", "signif")
  t_encoremieux <- t_mieux %>% 
    mutate(Type = c(rep("Direct", 14), rep("Indirect", 12),
                    rep("Total",14), rep("Mediator", 9))) %>%
    mutate_at(c("Effect", "Bias", "Std_Err", "Lower_CI", "Upper_CI"),
              as.numeric) %>%
    mutate_at(c("Type", "Var", "signif"), as.character) %>%
    mutate(Var = gsub(" ", "", Var))%>%
    mutate(Var = factor(Var, levels = c(
      "Lat", "Area", "Dist", "SLMP", "Elev", "GDP", "pop", "connect" ,"modif.change", 
      "static.modif", "CP", "SR.nat", "PD.nat","SES.PD.nat","FD.nat", "SES.FD.nat", 
      "SR.alien"), ordered = T))%>%
    mutate(signif = factor(signif, levels = c("*", ""), ordered = T))
  return(t_encoremieux)
}


fd_fig <- shape_tbl(t_eff_fd) %>% mutate(response = "FD")
pd_fig <- shape_tbl(t_eff_pd) %>% mutate(response = "PD")
cp_fd_fig <- shape_tbl_CP(t_eff_fd_cp) %>% mutate(response = "FD_CP")
cp_pd_fig <- shape_tbl_CP(t_eff_pd_cp) %>% mutate(response = "PD_CP")


eff_pd$Summary$SR.alien
eff_pd$Summary$PD.nat

# combine with and without CP
# FD
fd_and_fdcp <- bind_rows(fd_fig, cp_fd_fig) %>%
  filter(Type %in% c("Direct", "Indirect", "Total"))
# PD
pd_and_pdcp <- bind_rows(pd_fig, cp_pd_fig) %>%
  filter(Type %in% c("Direct", "Indirect", "Total"))
                      
# position of points
pd = position_dodge(.5)

# make the plots

p1bis <- ggplot(data = fd_and_fdcp, aes(x = Var, y = Effect, color = response,
                                       alpha = signif))+
  geom_abline(slope = 0, intercept = 0, linetype = 2)+
  geom_errorbar(aes(ymin=Lower_CI, ymax=Upper_CI),
                width=.4, position = pd, size = .6) +
  geom_point(position = pd, size=3, shape = 18) +
  scale_alpha_manual(values = c(1, .3))+
  scale_color_manual(values=c('turquoise', 'turquoise4'))+
  facet_wrap(~Type)+
  coord_flip()+
  theme_light()+
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12))
p2bis <- ggplot(data = pd_and_pdcp, aes(x = Var, y = Effect, color = response, 
                                    alpha = signif))+
  geom_abline(slope = 0, intercept = 0, linetype = 2)+
  geom_errorbar(aes(ymin=Lower_CI, ymax=Upper_CI),
                width=.4, position = pd, size = .6) +
  geom_point(position = pd, size=3, shape = 18) +
  scale_alpha_manual(values = c(1, .3))+
  scale_color_manual(values=c('springgreen', 'springgreen4'))+
  facet_wrap(~Type)+
  coord_flip()+
  theme_light()+
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12))
p1bis
p2bis

pdf

pdf("Fig/11_Effects_dir_indir_tot_fd_CP_noCP.pdf", 8, 5)
p1bis
dev.off()
pdf("Fig/11_Effects_dir_indir_tot_pd_CP_noCP.pdf", 8, 5)
p2bis
dev.off()

