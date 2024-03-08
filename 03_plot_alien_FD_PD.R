# Figure 1 with spatial vizualisation of data 
# adding PD and FD
rm(list=ls())


library(tidyverse)
library(sf)

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
  select(-c(varP, varT, HFI, GDP_mean, GDP_max, GDP_0_free, intact,
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

# scale PD & FD between 0 and 1 

all_scaled <- all %>% 
  mutate(PD_log_scale = (log(PD_alien)-min(log(PD_alien)))/
           (max(log(PD_alien))-min(log(PD_alien))),
         FD_log_scale = (log(fric_alien)-min(log(fric_alien)))/
           (max(log(fric_alien))-min(log(fric_alien)))) %>%
  mutate(PD_scale = (PD_alien-min(PD_alien))/
           (max(PD_alien)-min(PD_alien)),
         FD_scale = (fric_alien-min(fric_alien))/
           (max(fric_alien)-min(fric_alien)))

##################### Maps #####################

countries <- rnaturalearth::ne_countries(returnclass = "sf")

# with log_scaled values
logplot <- ggplot(countries)+
  geom_sf(color = NA, fill = "grey60", alpha = .7) + 
  # geom_point(data = all_scaled, 
  #            aes(x = Long, y = Lat, color = FD_log_scale, fill = PD_log_scale),
  #            shape = 21, stroke=2, size = 2)+
  geom_point(data = all_scaled, 
             aes(x = Long, y = Lat, color = FD_log_scale, size = PD_log_scale),
             shape = 1, stroke=1.3)+
  scale_colour_viridis_c() +
  xlab("") + ylab("")+
  theme_bw()+
  theme(legend.direction = "vertical", legend.position = "right")
logplot


# keep only realms with more than 10 islands for boxplots

realms_above10 <- sort(unique(all$Realm))[table(all$Realm) > 10]

all_long <- all_scaled %>%
  dplyr::select(ID, Realm, FD_scale, PD_scale) %>%
  pivot_longer(cols  = FD_scale:PD_scale,
               names_to = "metric",
               values_to = "value") %>% 
  dplyr::filter(Realm %in% realms_above10)

fd_pd <- ggplot(all_long, aes(Realm, value, color = Realm)) +
  geom_boxplot(size=.8, fill = NA, color = "grey50", outlier.color = NA) + 
  geom_jitter(height = 0, width = 0.2, alpha = 0.6, size = 2) + 
  scale_color_brewer(palette = "Paired")+
  xlab("")+ ylab("") +
  theme_classic(base_size = 15)+
  coord_flip()+
  facet_wrap(~metric, ncol = 2)+
  theme(legend.position = "none", strip.background = element_blank())

fd_pd

cor_pd_fd <- ggplot(all_scaled, aes(FD_scale, PD_scale)) +       
  geom_point(size=4, alpha=0.5, color = "grey50") +
  geom_smooth(method='lm', show.legend = T, size = 2, color = "grey10", se = T)+
  theme_classic(base_size = 15)+
  ylab("PD per island")+
  xlab("FD per island")

cor_pd_fd

library(patchwork)
#define the layout you want
layout <- "
AAAAAA
AAAAAA
AAAAAA
BBBBBB
BBBBBB
"

pdf("Fig/17_Fig_1_Map_FD_PD_world.pdf", 9,6)
logplot + fd_pd  + plot_layout(design = layout)
dev.off()

pdf("Fig/17_Fig_1_FD_PD_values.pdf", 7,3)
fd_pd
dev.off()


################ BASIC ANALYSES RESPONSE VARIABLES #################


summary(all$SES_fric_alien)
sum(abs(all$SES_fric_alien)>1.96)
sum((all$SES_fric_alien)>1.96)

summary(all$SES_PD_alien)
sum(abs(all$SES_PD_alien)>1.96)
sum((all$SES_PD_alien)>1.96)

cor(all$fric_alien, all$PD_alien)
cor.test(all$fric_alien, all$PD_alien)

cor.test(all_scaled$fric_alien, all_scaled$PD_alien)

plot(all$fric_alien, all$PD_alien)

plot(all_scaled$PD_log_scale, all_scaled$FD_log_scale)
cor.test(all_scaled$PD_log_scale, all_scaled$FD_log_scale)

colnames(all_scaled)

ggplot(all_scaled, aes(fric_alien, PD_alien))+
  geom_point()+
  geom_smooth()

ggplot(all_scaled, aes(FD_scale, PD_scale))+
  geom_point()+
  geom_smooth()

ggplot(all_scaled, aes(PD_log_scale, FD_log_scale))+
  geom_point()+
  geom_smooth()

ggplot(all_scaled, aes(SR_alien, PD_log_scale))+
  geom_point()+
  geom_smooth()
ggplot(all_scaled, aes(SR_alien, FD_log_scale))+
  geom_point()+
  geom_smooth()
ggplot(all_scaled, aes(SR_alien, PD_scale))+
  geom_point()+
  geom_smooth()
ggplot(all_scaled, aes(SR_alien, FD_scale))+
  geom_point()+
  geom_smooth()


cor.test(all_scaled$SR_alien, all_scaled$PD_alien)
cor.test(all_scaled$SR_alien, all_scaled$fric_alien)


dat <- read.csv2("Z:/THESE/6_Projects/biogeo_FD_alien/biogeo_alien_FD_PD/Data/Marino_et_al_DATA_407_ISL.csv")
colnames(dat)
ggplot(dat, aes(SR_alien, FD_alien))+
  geom_point()+
  geom_smooth()
ggplot(dat, aes(SR_alien, PD_alien))+
  geom_point()+
  geom_smooth()
ggplot(dat, aes(FD_alien, PD_alien))+
  geom_point()+
  geom_smooth()
