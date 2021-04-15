## you need to run import_data.R before running this script
source("scripts/MCMC_pg.R")

# these will not be changed
interaction_bool <- FALSE
neigh_bool <- TRUE
noNH_WHITE_bool <- TRUE
newseed_bool <- FALSE

burnin <- 500
thin <- 10
mcmc_niter <- 500*thin+burnin # with two chains we should have a total of 1000

############################### All data & noSMOTE ###############################
smote_bool <- FALSE

############### PRETERM ###############
outcome <- "PRETERM"

noRE_bool <- FALSE
prior <- "CAR"
source("scripts/philly_alldata.R")


noRE_bool <- FALSE
prior <- "independent"
source("scripts/philly_alldata.R")


noRE_bool <- TRUE
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
