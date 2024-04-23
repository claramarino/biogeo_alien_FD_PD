# Data

## Description of the main data file
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
01_islands_from_weigelt_et_al.rds: database of global islands from Weigelt et al. (2014)
07_oce_isl.rds: dataset containing all the contextual variable we retrieved for the oceanic islands contained in the Weigelt dataset
11_pd_+_ses_natives_world.rds

