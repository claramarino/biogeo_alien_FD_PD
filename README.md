

# Biogeography of alien birds
<!-- badges: start -->
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://choosealicense.com/licenses/mit/)
<!-- badges: end -->

Scripts and data for reproducing the results obtained by Marino et al. in the paper **The Anthropocene biogeography of alien birds on islands: drivers of their functional and phylogenetic diversities**.

## Description of the Data
All information related to data are detailed in `Data/README.md`. The main data file `Marino_et_al_DATA_407_ISL.csv` contains the name, location, and characteristics of all the 407 islands studied in the associated paper. Lines refer to islands.

## How to use the scripts
The script `00_create_clean_database_final.R` creates the cleaned data to use on the SEM analyses. Initial databases were too heavy for being loaded in the repository. Yet, we still provide the scripts that permit to calculate the functional an phylogenetic diversity (FD and PD) as well as the standardized effect sizes (SES-FD and SES-PD): `01_calculate_FD_SES_FD.R`, `02_calculate_PD_SES_PD.R`. The figures associated with spatial variation of, and relationships between, FD, PD, SES-FD, and SES-PD are coded in the script `03_plot_alien_FD_PD_and_SES.R`.

The scripts numbered from 10 to 13 perform the SEMs. Main SEMs presented in the paper are computed in the script `10_SEM_main_407_isl.R` while all other SEMs presented in Supplementary Material for robustness analyses are coded in scripts `11_SEM_subset_CP_96_isl.R`, `12_SEM_subset_CP_TSFI_83_isl.R`, `13_SEM_main_407_isl_no_human_var.R`.

Finally, the scripts 20, 21, and 22 plot the results of the SEM as direct, indirect and total effects.

## Citation

Please use the following citation:

> Marino C & Journiac L (2024) Data and codes for The Anthropocene biogeography of alien birds on islands: drivers of their functional and phylogenetic diversities. URL: https://github.com/claramarino/biogeo_alien_FD_PD.
