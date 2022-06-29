# geospatial-adverse-pregnancy

This repository contains R code to reproduce the results in the paper "Uncovering Patterns for Adverse Pregnancy Outcomes with Spatial Analysis: Evidence from Philadelphia" (Balocchi, Bai, Liu, Canelón, George, Chen and Boland, 2021+).

Please follow the following instructions to run the code.

1. Install the following R packages used for the data analysis and replication of figures and tables:
```
install.packages(c("mvtnorm", "BayesLogit", "sp", "maptools", "ggplot2", 
                    "ggmap", "dplyr", "rgeos", "pROC", "spdep", "bde", "xtable", "LaplacesDemon", "performanceEstimation"))
```

2. Import and merge the data patient level data and Philadelphia neighborhood data using `script/import_data.R`. Synthetic patient level data are available in `data/`, together with all the Philadelphia neighborhood data used for our analysis, merged into one file. Further details on how to download and merge the original data sources are provided in `data/`.

3. Use `script/meta_script.R` to run all or some of the data analyses considered in the paper. This corresponds to running the CAR model for both outcomes. We include also scripts for running alternative models, that were considered in the first version of the paper (for different prior distributions on the random effects, using the entire dataset or leave-one-out (LOO) by year). Further details are provided in `scripts/readme.md`. Output files from the analyses will be saved in `results/`. Additionally, `script/meta_script_supplement.R` provides the code to run the models reported in the supplementary material.

4. The scripts `create_figures.R`, `create_maps.R`, `create_tables.R` and `neigh_cluster_analysis.R` are used to reproduce the figures, maps and tables in the paper.
