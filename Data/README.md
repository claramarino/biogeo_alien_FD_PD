# Data
The data folder contains the main data file used for performing all the analyses from the associated paper (Marino_et_al_DATA_407_ISL.csv). However, we also provide intermediate data files that can be useful for understanding the process of final data obtention, or for exploring deeper the different robustness analyses we conducted in the Supplementary Material.

## Description of the main data file: Marino_et_al_DATA_407_ISL.csv
Data file Marino_et_al_DATA_407_ISL.csv contains the name, location, and characteristics of all the 407 islands studied in the associated paper. Lines refer to islands.

Column description:

1. Island nomenclature
- ID (character): island identifier from Weigelt et al. (2013)
- ARCHIP (character): archipelago name
- ISLAND (character): island name
- LONG (numeric): decimal longitude
- LAT (numeric): decimal latitude
- Realm (character): bioregion name

2. Geographic context (note that all variables are scaled and log-transformed if precised in the Material and Methods)
- Area (numeric): island area
- Dist (numeric): distance to the closest mainland
- SLMP (numeric): surrounding land mass
- Elev (numeric): island elevation
- Lat (numeric): island latitude

3. Human context (note that all variables are scaled and log-transformed if precised in the Material and Methods)
- pop (numeric): total population
- static_modif (numeric): Static habitat modification
- modif_change (numeric): Change in habitat modification
- GDP (numeric): gross domestic product
- connect (numeric): human connectivity 

4. Biotic context (note that all variables are scaled and log-transformed if precised in the Material and Methods)
- SR_nat (numeric): native taxonomic richness
- FD_nat (numeric): native functional richness
- SES_FD_nat (numeric): standardized effect size for native functional richness
- PD_nat (numeric): native phylogenetic richness
- SES_PD_nat (numeric): standardized effect size for native phylogenetic richness

5. Response variables (note that all variables are scaled and log-transformed if precised in the Material and Methods)
- SR_alien (numeric): alien taxonomic richness
- FD_alien (numeric): alien functional richness
- SES_FD_alien (numeric): standardized effect size for alien functional richness
- PD_alien (numeric): alien phylogenetic richness
- SES_PD_alien (numeric): standardized effect size for alien phylogenetic richness

## Description of other provided data files

### Island data (geographic and human contexts)
- 01_islands_from_weigelt_et_al.rds: database of global islands from Weigelt et al. (2014)
- 07_oce_isl.rds: dataset containing all the contextual variable we retrieved for the oceanic islands. All data sources are detailed in the Methods section of the associated paper and are also listed in the References below.

### PD and SES-PD of native and alien birds
- 11_pd_+_ses_natives_world.rds: output of the scripts that compute PD and SES-PD, for native birds
- 11_pd_by_isl_world.rds: output of the scripts that compute PD and SES-PD for alien birds, using as global pool for null models all birds in the world (n=10,862)
- 11_pd_by_isl_world_alien_only.rds: output of the scripts that compute PD and SES-PD for alien birds, using as global pool for null models alien birds only (n=952)

### FD and SES-FD of native and alien birds
- 13_fric_+_ses_natives_isl_over3.rds: output of the scripts that compute FD and SES-FD, for native birds
- 13_func_div_ses_world.rds: output of the scripts that compute FD and SES-FD for alien birds, using as global pool for null models all birds in the world (n=10,862)
- 13_func_div_ses_world_alien_only.rds: output of the scripts that compute FD and SES-FD for alien birds, using as global pool for null models alien birds only (n=952)


## References
Birdlife International & Handbook of the Birds of the World (2020). Bird species distribution maps of the world.

CIENSIN (2018). Gridded Population of the World, Version 4 (GPWv4): Population Count, Revision 11.

Dyer, E.E., Redding, D.W. & Blackburn, T.M. (2017b). The global avian invasions atlas, a database of alien bird distributions worldwide. Sci. Data, 4, 1–12.

IUCN. (2022). The IUCN Red List of Threatened Species.

Theobald, D.M., Kennedy, C., Chen, B., Oakleaf, J., Baruch-Mordo, S. & Kiesecker, J. (2020). Earth transformed: detailed mapping of global human modification from 1990 to 2017. Earth Syst. Sci. Data, 12, 1953–1972.

Tobias, J.A., Sheard, C., Pigot, A.L., Devenish, A.J.M., Yang, J., Neate-Clegg, M.H.C., et al. (2021). AVONET: morphological, ecological and geographical data for all birds. Ecol. Lett., 1–17.

Wang, T. & Sun, F. (2022). Global gridded GDP data set consistent with the shared socioeconomic pathways. Sci. Data, 9, 221.

Weigelt, P., Jetz, W. & Kreft, H. (2013). Bioclimatic and physical characterization of the world’s islands. Proc. Natl. Acad. Sci., 110, 15307–16342.

Other links:
- https://data.humdata.org/dataset/world-port-index
- http://www.partow.net/miscellaneous/airportdatabase/
