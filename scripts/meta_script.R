## you need to run import_data.R before running this script
library(dplyr)
library(pROC)
library(sp)
library(DMwR) # Needed for the SMOTE method

source("scripts/MCMC_pg.R")

# these will not be changed
interaction_bool <- FALSE		# do not consider interactions between patient-level and neighborhood level races variables
neigh_bool <- TRUE				# include neighborhood-level data among the covariates
noNH_WHITE_bool <- TRUE			# do not include the patient-level race variable for `White`, to avoid collinearity between covariates
newseed_bool <- FALSE			# consider a different seed to avoid running SMOTE with some neighborhoods having no data

burnin <- 500
thin <- 10
mcmc_niter <- 500*thin+burnin 	# with two chains we should have a total of 1000

############################### All data & noSMOTE ###############################
smote_bool <- FALSE

############### PRETERM ###############
outcome <- "PRETERM"

noRE_bool <- FALSE				# use a model with random effects
prior <- "CAR"
source("scripts/philly_alldata.R")


noRE_bool <- FALSE				# use a model with random effects
prior <- "independent"
source("scripts/philly_alldata.R")


noRE_bool <- TRUE				# use a model without random effects
source("scripts/philly_alldata.R")

############### STILLBIRTH ############### 
outcome <- "STILLBIRTH"

noRE_bool <- FALSE
prior <- "CAR"
source("scripts/philly_alldata.R")


noRE_bool <- FALSE
prior <- "independent"
source("scripts/philly_alldata.R")


noRE_bool <- TRUE
source("scripts/philly_alldata.R")

############################### All data & SMOTE (Only CAR now) ###############################

smote_bool <- TRUE

############### PRETERM ###############
outcome <- "PRETERM"

noRE_bool <- FALSE
prior <- "CAR"
source("scripts/philly_alldata.R")

############### STILLBIRTH ############### 
outcome <- "STILLBIRTH"
newseed_bool <- TRUE

noRE_bool <- FALSE
prior <- "CAR"
source("scripts/philly_alldata.R")
newseed_bool <- FALSE


############################### LOO & SMOTE (Only CAR now) ###############################
smote_bool <- TRUE

############### PRETERM ###############
outcome <- "PRETERM"

noRE_bool <- FALSE
prior <- "CAR"
source("scripts/philly_LOO.R")



############### STILLBIRTH ############### 
outcome <- "STILLBIRTH"
newseed_bool <- TRUE

noRE_bool <- FALSE
prior <- "CAR"
source("scripts/philly_LOO.R")
newseed_bool <- FALSE

############################### LOO & noSMOTE (Only CAR now) ###############################
smote_bool <- FALSE

############### PRETERM ###############
outcome <- "PRETERM"

noRE_bool <- FALSE
prior <- "CAR"
source("scripts/philly_LOO.R")



############### STILLBIRTH ############### 
outcome <- "STILLBIRTH"

noRE_bool <- FALSE
prior <- "CAR"
source("scripts/philly_LOO.R")
