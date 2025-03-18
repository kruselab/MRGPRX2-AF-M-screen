# MRGPRX2 virtual nanobody screen using AlphaFold-Multimer
This repository contains code used to perform the virtual screen described in [_In silico_ discovery of nanobody binders to a G-protein coupled receptor using AlphaFold-Multimer (Harvey _et al._ 2025)](https://www.biorxiv.org/content/10.1101/2025.03.05.640882v1).

### Installation guide
To install, simply clone repository or download files and source the files `scoring.R` and `evaluation.R`. Install time on a standard desktop computer is a few seconds.

### Instructions for use
Outputs from the ColabFold implementation of AlphaFold-Multimer (version 3) can be parsed and scored using the `run_pipeline()` command found in `scoring.R`. Functions related to the assessment of individual features' predictive ability can be found in `evaluation.R`. Details on the parameters used to generate ColabFold outputs can be found in the main text and methods of our preprint.

### Demo
A small example dataset is provided in `data/`. To score this dataset, run the command `run_pipeline("data/")`. Demo run time on a standard desktop computer is less than one minute.

Expected outputs:  
- A folder `res_pair_data/` containing CSV files, one per AF-M model, each of which contain information about interchain residue pairs in the specified model  
- A folder `score_data/` containing two files: `complex_scores.csv`, a summary file containing score information for each complex (1 row per AF-M run); and `model_scores.csv`, a summary file containing score information for each AF-M model (5 rows per AF-M run)

### System requirements
This code requires R as well as the following packages: `tidyverse`, `cowplot`, `RColorBrewer`, `jsonlite`, `micropan`, `Biostrings`, `GGally`, `plotROC`, `here`. This code has been tested using R version 4.2.2 on Mac OS Sequoia version 15.3.2.
