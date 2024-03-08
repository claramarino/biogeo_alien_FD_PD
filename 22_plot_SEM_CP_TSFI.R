# SEM figure & tables
# for models with CP and TSFI

rm(list=ls())

library(semEff)
library(tidyverse)


# load models

# FD without CP/TSFI
mod_fd <- readRDS("Output/11_SEM_noCP_FD_model.rds")
eff_fd <- readRDS("Output/11_SEM_noCP_FD_boot_eff.rds")
# PD without CP/TSFI
mod_pd <- readRDS("Output/11_SEM_noCP_PD_model.rds")
eff_pd <- readRDS("Output/11_SEM_noCP_PD_boot_eff.rds")

# FD with CP t& TFI
mod_fd_CT <- readRDS("Output/12_SEM_CP_TSFI_FD_model.rds")
eff_fd_CT <- readRDS("Output/12_SEM_CP_TSFI_FD_boot_eff.rds")
# PD with CP & TFI
mod_pd_CT <- readRDS("Output/12_SEM_CP_TSFI_PD_model.rds")
eff_pd_CT <- readRDS("Output/12_SEM_CP_TSFI_PD_boot_eff.rds")

# Extract effects as data.frame for figures

t_eff_fd <- as.data.frame(eff_fd$Summary$FD.alien)
t_eff_pd <- as.data.frame(eff_pd$Summary$PD.alien)
t_eff_fd_ct <- as.data.frame(eff_fd_CT$Summary$FD.alien)
t_eff_pd_ct <- as.data.frame(eff_pd_CT$Summary$PD.alien)

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
      "Lat", "Area", "Dist", "SLMP", "Elev", "GDP", "pop", "TSFI", "connect" ,"modif.change", 
      "static.modif", "CP","SR.nat", "PD.nat","SES.PD.nat","FD.nat", "SES.FD.nat", 
      "SR.alien"), ordered = T))%>%
    mutate(signif = factor(signif, levels = c("*", ""), ordered = T))
  return(t_encoremieux)

}

# use a function to reshape tables with CP AND TFI
shape_tbl_CT <- function(t_pourri){
  t_mieux <- t_pourri[-c(1,17,31,47), -c(3,5,7,9,12)]
  
  names(t_mieux) <- c("Type", "Var", "Effect", "Bias", "Std_Err", "Lower_CI", "Upper_CI", "signif")
  t_encoremieux <- t_mieux %>% 
    mutate(Type = c(rep("Direct", 15), rep("Indirect", 13),
                    rep("Total",15), rep("Mediator", 10))) %>%
    mutate_at(c("Effect", "Bias", "Std_Err", "Lower_CI", "Upper_CI"),
              as.numeric) %>%
    mutate_at(c("Type", "Var", "signif"), as.character) %>%
    mutate(Var = gsub(" ", "", Var))%>%
    mutate(Var = factor(Var, levels = c(
      "Lat", "Area", "Dist", "SLMP", "Elev", "GDP", "pop", "TSFI", "connect" ,"modif.change", 
      "static.modif", "CP", "SR.nat", "PD.nat","SES.PD.nat","FD.nat", "SES.FD.nat", 
      "SR.alien"), ordered = T))%>%
    mutate(signif = factor(signif, levels = c("*", ""), ordered = T))
  return(t_encoremieux)
}

fd_fig <- shape_tbl(t_eff_fd) %>% mutate(response = "FD")
pd_fig <- shape_tbl(t_eff_pd) %>% mutate(response = "PD")
ct_fd_fig <- shape_tbl_CT(t_eff_fd_ct) %>% mutate(response = "FD_CT")
ct_pd_fig <- shape_tbl_CT(t_eff_pd_ct) %>% mutate(response = "PD_CT")


eff_pd$Summary$SR.alien
eff_pd$Summary$PD.nat

# position of points
pd = position_dodge(.5)

# combine with and without CP/TSFI
# FD
fd_and_fdct <- bind_rows(fd_fig, ct_fd_fig) %>%
  filter(Type %in% c("Direct", "Indirect", "Total"))
# PD
pd_and_pdct <- bind_rows(pd_fig, ct_pd_fig) %>%
  filter(Type %in% c("Direct", "Indirect", "Total"))

# make the plots

p1bis <- ggplot(data = fd_and_fdct, aes(x = Var, y = Effect, color = response,
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

p2bis <- ggplot(data = pd_and_pdct, aes(x = Var, y = Effect, color = response, 
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



# save plots as pdf

pdf("Fig/12_Effects_dir_indir_tot_fd_CP_TSFI.pdf", 8, 5)
p1bis
dev.off()
pdf("Fig/12_Effects_dir_indir_tot_pd_CP_TSFI.pdf", 8, 5)
p2bis
dev.off()


