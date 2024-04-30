#savefiles 
saveRDS(all_scaled, "Data/data1_maps_with_ggplot.RDS")

# load datasets
data1 <- readRDS("Data/data1_maps_with_ggplot.RDS")

library(tidyverse) # contains ggplot
# for handling spatial objects: package "sf"
# install.packages("sf")

# *short introduction to spatial packages in R*

