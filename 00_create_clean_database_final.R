# create clean public data
# fd, pd, ses, insular contexts (for SEM)
rm(list=ls())


library(tidyverse)

################ Load data ###############

# contextual variables
var <- readRDS("Data/07_oce_isl.rds") 

# alien diversities
fric_alien <- readRDS("Data/13_func_div_ses_world_alien_only.rds")
pd_alien <- readRDS("Data/11_pd_by_isl_world_alien_only.rds")

# native diversities
fric_nat <- readRDS("Data/13_fric_+_ses_natives_isl_over3.rds")
pd_nat <- readRDS("Data/11_pd_+_ses_natives_world.rds")

# islands from Weigelt et al 2013
wi_ok <- readRDS("Data/01_islands_from_weigelt_et_al.rds")

# time of first introduction +
# colonization pressure info for a subset of islands
tsi <- readRDS("Data/R1_TSI_83_isl.rds")
cp <- readRDS("Data/R1_Col_pressure_96_isl.rds")


######### Clean and join tables ############


var_ok <- left_join(var %>% 
  mutate(connect = nb_ap + nb_ap_buff_100) %>% # transform connect as slmp 
  # remove var not used in models
  select(-c(varP, varT, HFI, GDP_mean, GDP_max, GDP_0_free, intact, Long,
            GMMC, nb_ap, nb_ap_buff_50, nb_ap_buff_100, dist_ap)) %>%
  rename(SR_nat = native_sp_rich,
         GDP = GDP_sum) %>% 
  mutate(ID = as.character(ID)),
  left_join(
    cp %>% 
      mutate(ID = as.character(ID)) %>%
      rename(CP=n),
    tsi %>% 
      mutate(ID = as.character(ID),
             TSFI = 2020 - TSI)) %>%
    select(-TSI))


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
              "PD_nat", "SR_alien", "SR_nat", "CP", "TSFI"), log) %>%
  mutate_at(c("pop", "Elev","GDP", "connect"), ~log(.+1)) %>%
  mutate(Lat = abs(Lat)) %>%
  # scale variables before SEM
  mutate_if(is.numeric, ~as.numeric(scale(.)))


table(all_ok$Realm)

colSums(is.na(all_ok))



all_isl <- inner_join(wi_ok, all_ok) %>%
  rename(FD_alien = fric_alien, 
         SES_FD_alien = SES_fric_alien,
         FD_nat = fric_nat,
         SES_FD_nat = SES_fric_nat) %>%
  select(ID, ARCHIP, ISLAND, LONG, LAT, Realm,
         Area, Dist, SLMP, Elev, Lat,
         pop, static_modif, modif_change, GDP, connect, CP, TSFI,
         SR_nat, FD_nat, SES_FD_nat, PD_nat, SES_PD_nat,
         SR_alien, FD_alien, SES_FD_alien, PD_alien, SES_PD_alien)


colnames(all_isl)

write.csv2(all_isl, "Data/Marino_et_al_DATA_407_ISL.csv", row.names = F)




#### Check the influence of pool for null models (SES) ####

all_isl <- read.csv("Data/Marino_et_al_DATA_407_ISL.csv", sep = ";")

# pool 1 => null models are drawn within all birds (n=10,862)
fd_alien_pool1 <- readRDS("Data/13_func_div_ses_world.rds") %>%
  .[.$isl%in% all_isl$ID, ]
pd_alien_pool1 <- readRDS("Data/11_pd_by_isl_world.rds") %>%
  .[.$isl%in% all_isl$ID, ]

# pool 2 => null models are drawn within alien birds only (n=952)
pd_alien_pool2 <- readRDS("Data/11_pd_by_isl_world_alien_only.rds") %>%
  .[.$isl%in% all_isl$ID, ]
fd_alien_pool2 <- readRDS("Data/13_func_div_ses_world_alien_only.rds")

# Compare the two pools
pd_all <- left_join(pd_alien_pool1 %>% rename(ses_PD_all = ses_PD,
                                        PD_mean1 = PD_mean),
                    pd_alien_pool2 %>% rename(ses_PD_alien = ses_PD,
                                           PD_mean2 = PD_mean))

fd_all <- left_join(fd_alien_pool1 %>% rename(ses_FD_all = SES_fric,
                                              FD_mean1 = fric),
                    fd_alien_pool2 %>% rename(ses_FD_alien = SES_fric,
                                              FD_mean2 = fric))

ggplot(pd_all, aes(x=ses_PD_all, y = ses_PD_alien))+
  geom_point()
hist(pd_all$ses_PD_all)
hist(pd_all$ses_PD_alien)
t.test(pd_all$ses_PD_alien)
# t = -13.562, df = 406, p-value < 2.2e-16, mean = -0.8937138
cor.test(pd_all$ses_PD_all, pd_all$ses_PD_alien)
# cor = 0.9784464 

ggplot(fd_all, aes(x=ses_FD_all, y = ses_FD_alien))+
  geom_point()
hist(fd_all$ses_FD_all)
hist(fd_all$ses_FD_alien)
t.test(fd_all$ses_FD_alien)
# t = -33.579, df = 406, p-value < 2.2e-16, mean = -0.9188079
cor.test(fd_all$ses_FD_all, fd_all$ses_FD_alien)
# cor = 0.9556737 


# the two pools give very similar values for final SES FD/PD
# we select the more conservative, ie the pool 2 with all alien birds