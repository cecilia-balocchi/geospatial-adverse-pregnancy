## Importing the data

The script `import_data.R` loads the fake pregnancies dataset, merges it with the Philadelphia neighborhood data and performs the needed transformations to get the dataset ready for running the models.

## Running the models

The script `meta_script.R` runs all the models considered. It imports code from other scripts:
- `MCMC_pg.R` contains the function to run the MCMC using Polya-Gamma (pg),
- `philly_train_test.R` runs the MCMC with testing set defined by the last year of data (default)
- `philly_alldata.R` runs the MCMC when all the data is used (no leave one out by year) (optional, discarded)
- `philly_LOO.R` runs the MCMC when the leave one out by year is used. (optional, discarded)

Several options are set before running each model, which are used to select the right code within `philly_train_test.R`:
- `outcome`, which can be "PRETERM" or "STILLBIRTH"
- `noRE_bool`, which is set to FALSE when we want to run a model with random effects
- `prior`, which can be "CAR" or "independent", selects the prior for the random effects and is used when `noRE_bool` is equal to FALSE.

Similarly, the script `meta_script_supplement.R` runs the bayesian models which are reported in the supplement. 

## Creating tables and figures

The scripts `create_figures.R`, `create_maps.R` and `create_tables.R` creates the figures, maps and tables reported in the manuscript and supplement, except for those related to the neighborhood cluster analysis, which are created in `neigh_cluster_analysis.R`.
The script `functions_scale.R` was downloaded from http://egallic.fr/en/scale-bar-and-north-arrow-on-a-ggplot2-map/ and used to add a bar with a scale for maps.