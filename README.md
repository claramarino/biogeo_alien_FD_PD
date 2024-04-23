# Biogeography of alien birds
Scripts and data for reproducing the results obtained by Marino et al. in the paper "The Anthropocene biogeography of functional and phylogenetic diversities of alien birds"

## Description of the Data and file structure
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

## How to use the scripts
The script 00_create_clean_database_final.R creates the cleaned data to use on the SEM analyses. Initial databases were too heavy for being loaded in the repository. Yet, we still provide the scripts that permit to calculate the functional an phylogenetic diversity (FD and PD) as well as the standardized effect sizes (SES-FD and SES-PD): 01_calculate_FD_SES_FD.R, 02_calculate_PD_SES_PD.R. The figures associated with spatial varaition of and relationships between FD, PD, SES-FD, and SES-PD are coded in the script 03_plot_alien_FD_PD_and_SES.R.

The scripts numbered from 10 to 13 perform the SEMs. Main SEMs presented in the paper are computed in the script 10_SEM_main_407_isl.R while all other SEMs presented in Supplementary Material for robustness analyses are coded in scripts 11_SEM_subset_CP_96_isl.R, 12_SEM_subset_CP_TSFI_83_isl.R, 13_SEM_main_407_isl_no_human_var.R.

Finally, the scripts 20, 21, and 22 plot the results of the SEM as direct, indirect and total effects.
