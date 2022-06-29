## you need to run import_data.R before running this script
library(dplyr)
library(pROC)
library(sp)
library(LaplacesDemon)  # for MCSE and ESS

## packages that are not needed for replication, but that we used at some point
# library(MLmetrics) 				# to compute PR-AUC (uncomment line 494 in philly_train_test.R to compute PRAUC)

source("scripts/MCMC_pg.R")
source("scripts/functions_for_race.R")

interaction_bool <- FALSE		# do not consider interactions between patient-level and neighborhood level races variables
neigh_bool <- TRUE				# include neighborhood-level data among the covariates
newseed_bool <- FALSE			# consider a different seed to avoid running SMOTE with some neighborhoods having no data
smote_bool <- FALSE				# do not use SMOTE (discarded because our focus in not only on predictions)
reweight_bool <- FALSE			# do not use a reweighted likelihood (discarded because our focus in not only on predictions)
alpha_vs_beta_bool <- TRUE 		# determines which kind of reweight. Does not matter because reweight_bool is FALSE

burnin <- 500
thin <- 5
mcmc_niter <- 1000*thin+burnin	# with two chains we should have a total of 2000 (and approximately 1600 ESS)

##############################################################
#################### Models in main paper ####################
##############################################################

########################## CAR model #########################

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

# To replicate the results in the paper, you do not need to run scripts past this point.

##############################################################
############ Additional models of possible interest ##########
##############################################################

################## independent RE model ######################

# we compared the performance of CAR model with independent RE model

###### PRETERM ######
outcome <- "PRETERM"

noRE_bool <- FALSE				# use a model with random effects
prior <- "independent"			# use a conditionally iid prior on the random effects
source("scripts/philly_train_test.R")

###### STILLBIRTH ###### 
outcome <- "STILLBIRTH"

noRE_bool <- FALSE
prior <- "independent"
source("scripts/philly_train_test.R")

####################### no RE model ###########################

# we initially considered a model without Random Effects

###### PRETERM ######
outcome <- "PRETERM"

noRE_bool <- TRUE				# use a model without random effects
source("scripts/philly_train_test.R")

###### STILLBIRTH ###### 
outcome <- "STILLBIRTH"

noRE_bool <- TRUE
source("scripts/philly_train_test.R")

########################### LOO ###############################

# we considered a Leave-One-Out (LOO) approach where each year is iteratively used as testing set
# any prior and outcome can be used here

outcome <- "PRETERM"
noRE_bool <- FALSE
prior <- "CAR"
source("scripts/philly_LOO.R")
