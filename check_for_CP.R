# gavia data
library(tidyverse)

# open table
gavia = read.csv("Z:/THESE/5_Data/Alien_data/GAVIA_tables/GAVIA_main_data_table.csv")

colnames(gavia)
length(unique(gavia$SpeciesID))

#check nb of species with polygons: all species that have an established pop 
gavia_map <- gavia %>% # 921 species in total
  filter(RangeMap == "Mapped") %>% # 720 species were mapped
  filter(StatusCat=="Established") # 362 with a current established population


length(unique(gavia_map$SpeciesID))
# 362 with a current established population with Mapped info
length(unique(gavia %>% filter(StatusCat=="Established") %>% pull(SpeciesID)))
# 419 with a current established population (no shp for 57 species)


# what are those 57 species?
no_map <- gavia %>%
  filter(StatusCat=="Established") %>%
  filter(RangeMap != "Mapped") %>%
  filter(! SpeciesID %in% unique(gavia_map$SpeciesID))
length(unique(no_map$SpeciesID))

# introduced in how many places in average?
nb_loc <- no_map %>%
  group_by(SpeciesID) %>%
  count()

summary(nb_loc)
table(nb_loc$n)

table(no_map$LandType)

no_map %>% filter(LandType=="Oceanic.Island")  
  
  

table(gavia$Island, gavia$LandType)
gavia_oi <- gavia %>% 
  filter(LandType=="Oceanic.Island" | LandType=="" & Island=="TRUE") %>%
  filter(CPrecord=="TRUE")
length(unique(gavia_oi$SpeciesID))

sum(is.na(gavia_oi$IntroducedDateGrouped))
sum(is.na(gavia$IntroducedDateGrouped))


dat <- read.csv2("Z:/THESE/6_Projects/biogeo_FD_alien/biogeo_alien_FD_PD/Data/Marino_et_al_DATA_407_ISL.csv")


# isl wglt shp 
library(sf)
isl <- st_read("Z:/THESE/6_Projects/biogeo_FD_alien/Data/Islands_Weigelt_reparees.shp")
colnames(isl)
head(isl)

length(unique(gavia_oi$CountryName))
unique(gavia_oi$CountryName)
length(unique(gavia_oi$AreaName1))
length(unique(gavia_oi$AreaName2))


tbl <- read.csv("Z:/THESE/5_Data/Distribution_spatiale/Weigelt_isl_database/Weigelt_etal_2013_PNAS_islanddata.csv")

actri <- st_read("Z:/THESE/5_Data/Distribution_spatiale/Alien_species/GAVIA_rangemaps/GAVIA_Acridotheres_tristis.shp")
plot(actri)


####
# how many archipels ?
arch <- gavia_oi %>% 
  mutate(archip = paste(CountryName, AreaName1, sep = " "))
arch$archip <- tolower(arch$archip)

unique(arch$archip)

dat_count <- left_join(dat, tbl %>% select(ID, Country))

################
# go with manual matching of the databases

openxlsx::write.xlsx(dat_count %>% select(ID, Country, ARCHIP, ISLAND),
                     file = "Z:/THESE/6_Projects/biogeo_FD_alien/biogeo_alien_FD_PD/Output/407_isl_to_Match_gavia_CP.xlsx")

openxlsx::write.xlsx(gavia_oi %>% select(RecordID, CountryName:Realm),
                     file = "Z:/THESE/6_Projects/biogeo_FD_alien/biogeo_alien_FD_PD/Output/Gavia_CP_oceanic_isl.xlsx")

################
gav_to_match <- openxlsx::read.xlsx(paste0(
  "Z:/THESE/6_Projects/biogeo_FD_alien/biogeo_alien_FD_PD/Output/",
  "Gavia_CP_oceanic_isl_manual_edit.xlsx"
))

dat_to_match <- openxlsx::read.xlsx(paste0(
  "Z:/THESE/6_Projects/biogeo_FD_alien/biogeo_alien_FD_PD/Output/",
  "407_isl_to_Match_gavia_CP_manual_edit.xlsx"
))


table(dat_to_match$in_gavia)
length(unique(gav_to_match$ID))
length(unique(paste0(gav_to_match$ID, gav_to_match$Country, gav_to_match$ARCHIP)))
hist(dat$SR_alien)

isl_w_cp <- gav_to_match %>%
  group_by(ID) %>%
  count() %>%
  filter(!is.na(ID))
# 96 islands for which we have colonization pressure
hist(isl_w_cp$n, breaks = 15)

dat_cp <- left_join(dat, isl_w_cp)
plot(dat_cp$SR_alien, dat_cp$n)
ggplot(dat_cp, aes(x=SR_alien, y = n))+
  geom_point()+
  geom_smooth()

colnames(gavia)


tsi <- left_join(gav_to_match, gavia %>% select(RecordID, IntroducedDateGrouped))
min_tsi <- tsi %>%
  filter(!is.na(IntroducedDateGrouped)) %>%
  group_by(ID) %>%
  summarize(TSI = min(IntroducedDateGrouped)) %>%
  filter(!is.na(ID)) %>% filter(!is.na(TSI))

cp_tsi <- left_join(min_tsi, isl_w_cp)
cor.test(cp_tsi$n, cp_tsi$TSI)
hist(cp_tsi$TSI)


plot(cp_tsi$TSI, cp_tsi$n)

#save the colonization pressure info in a table for analyses
saveRDS(isl_w_cp, "Data/R1_Col_pressure_96_isl.rds")
# save tsi 
saveRDS(min_tsi, "Data/R1_TSI_83_isl.rds")



dat_cp_tsi <- left_join(dat, cp_tsi)
plot(dat_cp_tsi$TSI, dat_cp_tsi$SR_alien)

ggplot(dat_cp_tsi, aes(x=TSI, y = n))+
  geom_point() +
  theme_classic()+
  xlab("Date of first alien bird introduction")+
  ylab("Colonization pressure")

#saveRDS(dat_cp_tsi, "Data/R1_CP_96_TFI_80_isl.rds")


colnames(gavia)

dat_cp_tsi <- readRDS("Data/R1_CP_96_TFI_80_isl.rds")

#remove all islands with less than 3 CP because we must have at least 4 species on all islands



# aggregate by archip?

# remove all records not in our db
gav_clean <- gav_to_match %>%
  filter(!(is.na(ID) & is.na(Country) & is.na(ARCHIP)))

# for all records with an isl ID, attribute country and archipel
for(i in 1:nrow(gav_clean)){
  # i=12
  id = gav_clean$ID[i]
  country = gav_clean$Country[i]
  archip = gav_clean$ARCHIP[i]
  
  if(!is.na(id)){ 
    gav_clean$Country[i] <- dat_to_match$Country[dat_to_match$ID==id]
    gav_clean$ARCHIP[i] <- dat_to_match$ARCHIP[dat_to_match$ID==id]
  } else {
    if(!is.na(country) & is.na(archip)){
      if(length(unique(dat_to_match$ARCHIP[dat_to_match$Country==country]))==1){
        gav_clean$ARCHIP[i] <- unique(dat_to_match$ARCHIP[dat_to_match$Country==country])}
    }
  }
}

arch_w_cp <- gav_clean %>%
  mutate(count_arch = paste0(Country, ARCHIP)) %>%
  group_by(count_arch) %>%
  summarise(n_arch = n())

dat_cp_all <- left_join(left_join(dat_cp, dat_to_match %>% 
                                    select(ID, Country)) %>% 
                          mutate(count_arch = paste0(Country, ARCHIP)),
                        arch_w_cp)


summary(dat_cp_all$n_arch)
summary(dat_cp_all$n)

ggplot(data = dat_cp_all)+
  geom_point(aes(x=SR_alien, y = n_arch), color = "red")+
  geom_smooth(aes(x=SR_alien, y = n_arch), color = "red")+
  geom_point(aes(x=SR_alien, y = n), color = "blue")+
  geom_smooth(aes(x=SR_alien, y = n), color = "blue")

