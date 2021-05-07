# geospatial-adverse-pregnancy

This repository contains R code to reproduce the results in the paper "A Bayesian Hierarchical Modeling Framework for Geospatial Analysis of Adverse Pregnancy Outcomes" (Balocchi, Bai, Liu, Canelón, George, Chen and Boland, 2021+).

Please follow the following instructions to run the code.

1. Install the following R packages used for the data analysis and replication of figures and tables:
```
install.packages(c("mvtnorm", "BayesLogit", "sp", "maptools", "ggplot2", "ggmap", "dplyr", "rgeos", "pROC", "DMwR", "spdep"))
```

2. Import and merge the data patient level data and Philadelphia neighborhood data using `script/import_data.R`. Synthetic patient level data are available in `data/`, together with all the Philadelphia neighborhood data used for our analysis, merged into one file. Further details on how to download and merge the original data sources are provided in `data/`.

3. Use `script/meta_script.R` to run all or some of the data analyses considered in the paper, for different outcomes, different prior distributions on the random effects, using SMOTE or the original dataset, using the entire dataset or leave-one-out by year. Further details are provided in `scripts/readme.md`. Output files from the analyses will be saved in `results/`.

4. The scripts `create_figures.R`, `create_maps.R`, `create_tables.R` are used to reproduce the figures, maps and tables in the paper.
