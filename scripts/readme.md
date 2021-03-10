## Importing the data

The script `import_data.R` loads the fake pregnancies dataset, merges it with the philadelphia neighborhood data and performs the needed transformations to get the dataset ready for running the models.

## Running the models

The script `meta_script.R` runs all the models considered. It imports code from other scripts:
- `MCMC_pg_nogamma_rho.R` contains the function to run the MCMC,
- `philly_alldata.R` runs the MCMC when all the data is used (no leave one out by year),
- `philly_LOO.R` runs the MCMC when the leave one out by year is used.
Several options are set before running each model, which are used to select the right code within `philly_alldata.R` and `philly_LOO.R`:
- `smote_bool`, which is set to TRUE when we want to run a model using SMOTE
- `outcome`, which can be "PRETERM" or "STILLBIRTH"
- `noRE_bool`, which is set to FALSE when we want to run a model with random effects
- `prior`, which can be "CAR" or "independent", selects the prior for the random effects and is used when `noRE_bool` is equal to FALSE.
(interaction_bool, neigh_bool, noNH_WHITE_bool, newseed_bool)

## Creating tables and figures

The script `create_figures.R` creates the illustration for the bayes-p and the plot of the log-odds (Figures .. and .. in the manuscript).
The script `create_maps.R` creates the maps of the proportions of pregnancies events and of the predicted probabilities. (TODO complete create_maps.R with the map of proportions).

The script `create_tables.R` creates the table reported in the manuscript.
(TODO decide whether to keep or not the neighborhood comparison table)
