#savefiles 
saveRDS(all_scaled, "Data/data1_maps_with_ggplot.RDS")

# load datasets
data1 <- readRDS("Data/data1_maps_with_ggplot.RDS")

library(tidyverse) # contains ggplot
# for handling spatial objects: package "sf"
# install.packages("sf")

# *short introduction to spatial packages in R*

# get the background
# usefull package: rnaturalearth => obtain world backgrounds in sf format
# install.packages("rnaturalearth")

countries <- rnaturalearth::ne_countries(returnclass = "sf")
# you can adjust the scale
# or select only some countries or a continent
europe <- rnaturalearth::ne_countries(continent = "europe", returnclass = "sf")

# plot simple feature
ggplot(countries)+
  geom_sf()


ggplot(countries)+
  geom_sf(color = NA, # remove boarders
          fill = "grey60", 
          alpha = .7)+
  theme_bw()

# the sf object can contain a lot of attributes
str(countries)
# you can use them for mapping

# one color per continent
ggplot(countries)+
  geom_sf(
    # pass arguments into the aes() function
    aes(fill = continent),
    color = NA, # remove boarders
    #fill = "grey60", 
    alpha = .7)+
  theme_bw()

# color gradient for population density
ggplot(countries)+
  geom_sf(
    # pass arguments into the aes() function
    aes(fill = pop_rank),
    color = NA)+
  scale_fill_viridis_c()+
  theme_bw()



# add points from your study to the background
# the info we want to plot: 
# functional and phylogenetic diversity of alien birds on oceanic islands

# start with functional diversity

ggplot(countries)+
  # plot the background
  geom_sf(color = NA, fill = "grey60", alpha = .7) + 
  # add your points
  geom_point(data = all_scaled, # specify the dataset 
             aes(x = Long, 
                 y = Lat, 
                 color = FD_log_scale), # color will change according to FD
             shape = 19)+
  scale_colour_viridis_c() +
  xlab("") + ylab("")+
  theme_bw()+
  theme(legend.direction = "vertical", legend.position = "right")

ggplot(countries)+
  # plot the background
  geom_sf(color = NA, fill = "grey60", alpha = .7) + 
  # add your points
  geom_point(data = all_scaled, # specify the dataset 
             aes(x = Long, 
                 y = Lat, 
                 color = FD_log_scale, # color will change according to FD
                 size = FD_log_scale), # size will also change
             shape = 1, stroke = 1.3)+
  scale_colour_viridis_c() +
  xlab("") + ylab("")+
  theme_bw()+
  theme(legend.direction = "vertical", legend.position = "right")


# now add the phylogenetic diversity : replace the size parameter
ggplot(countries)+
  # plot the background
  geom_sf(color = NA, fill = "grey60", alpha = .7) + 
  # add your points
  geom_point(data = all_scaled, # specify the dataset 
             aes(x = Long, 
                 y = Lat, 
                 color = FD_log_scale, # color will change according to FD
                 size = PD_log_scale), # size will change according to PD
             shape = 1, stroke = 1.3)+
  scale_colour_viridis_c() +
  xlab("") + ylab("")+
  theme_bw()+
  theme(legend.direction = "vertical", legend.position = "right")

