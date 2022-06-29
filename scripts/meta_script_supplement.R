## you need to run import_data.R before running this script
library(dplyr)
library(pROC)
library(sp)
library(performanceEstimation) 	# Needed for the SMOTE method
library(coda)			# ?
library(LaplacesDemon)  # ?

source("scripts/MCMC_pg.R")
source("scripts/functions_for_race.R")

interaction_bool <- FALSE		# do not consider interactions between patient-level and neighborhood level races variables
neigh_bool <- TRUE				# include neighborhood-level data among the covariates
newseed_bool <- FALSE			# consider a different seed to avoid running SMOTE with some neighborhoods having no data
smote_bool <- FALSE				# this will be changed to TRUE before running models with SMOTE
reweight_bool <- FALSE			# this will be changed to TRUE before running models with reweighted likelihood
alpha_vs_beta_bool <- TRUE 		# determines which kind of reweight. TRUE corresponds to downweighting the controls and is faster (recommended). 

burnin <- 500
thin <- 5
mcmc_niter <- 1000*thin+burnin	# with two chains we should have a total of 2000 (and approximately 1600 ESS)

####################### CAR model & SMOTE ######################
smote_bool <- TRUE

###### PRETERM ######
outcome <- "PRETERM"

noRE_bool <- FALSE				# use a model with random effects
prior <- "CAR"					# use a CAR prior on the random effects
source("scripts/philly_train_test.R")

###### STILLBIRTH ###### 
outcome <- "STILLBIRTH"

noRE_bool <- FALSE
prior <- "CAR"
source("scripts/philly_train_test.R")

smote_bool <- FALSE
##################### CAR model & reweight ####################
reweight_bool <- TRUE

###### PRETERM ######
outcome <- "PRETERM"

noRE_bool <- FALSE				# use a model with random effects
prior <- "CAR"					# use a CAR prior on the random effects
source("scripts/philly_train_test.R")

###### STILLBIRTH ###### 
outcome <- "STILLBIRTH"

noRE_bool <- FALSE
prior <- "CAR"
source("scripts/philly_train_test.R")
