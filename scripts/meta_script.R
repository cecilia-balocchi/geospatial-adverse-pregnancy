fldr="/Users/bolandm/Box\ Sync/Geographic_analysis/"
source(paste(fldr,"scripts/MCMC_pg_nogamma_rho.R", sep=""))

# these will not be changed
interaction_bool <- FALSE
neigh_bool <- TRUE
noNH_WHITE_bool <- TRUE
newseed_bool <- FALSE

# test first
burnin <- 500
thin <- 10
mcmc_niter <- 500*thin+burnin # with two chains we should have a total of 1000

############################### All data & noSMOTE ###############################
smote_bool <- FALSE

############### PRETERM ###############
outcome <- "PRETERM" # "PRETERM","STILLBIRTH","CSECTION"

noRE_bool <- FALSE
prior <- "CAR" # "CAR" "independent"
source(paste(fldr,"scripts/philly_alldata.R", sep=""))


noRE_bool <- FALSE
prior <- "independent" # "CAR" "independent"
source(paste(fldr,"scripts/philly_alldata.R", sep=""))


noRE_bool <- TRUE
source(paste(fldr,"scripts/philly_alldata.R", sep=""))

############### STILLBIRTH ############### 
outcome <- "STILLBIRTH" # "PRETERM","STILLBIRTH","CSECTION"

noRE_bool <- FALSE
prior <- "CAR" # "CAR" "independent"
source(paste(fldr,"scripts/philly_alldata.R", sep=""))


noRE_bool <- FALSE
prior <- "independent" # "CAR" "independent"
source(paste(fldr,"scripts/philly_alldata.R", sep=""))


noRE_bool <- TRUE
source(paste(fldr,"scripts/philly_alldata.R", sep=""))

############################### All data & SMOTE (Only CAR now) ###############################

smote_bool <- TRUE

############### PRETERM ###############
outcome <- "PRETERM" # "PRETERM","STILLBIRTH","CSECTION"

noRE_bool <- FALSE
prior <- "CAR" # "CAR" "independent"
source(paste(fldr,"scripts/philly_alldata.R", sep=""))

############### STILLBIRTH ############### 
outcome <- "STILLBIRTH" # "PRETERM","STILLBIRTH","CSECTION"
newseed_bool <- TRUE

noRE_bool <- FALSE
prior <- "CAR" # "CAR" "independent"
source(paste(fldr,"scripts/philly_alldata.R", sep=""))
newseed_bool <- FALSE


############################### LOO & SMOTE (Only CAR now) ###############################
smote_bool <- TRUE

############### PRETERM ###############
outcome <- "PRETERM" # "PRETERM","STILLBIRTH","CSECTION"

noRE_bool <- FALSE
prior <- "CAR" # "CAR" "independent"
source(paste(fldr,"scripts/philly_LOO.R", sep=""))



############### STILLBIRTH ############### 
outcome <- "STILLBIRTH" # "PRETERM","STILLBIRTH","CSECTION"
newseed_bool <- TRUE

noRE_bool <- FALSE
prior <- "CAR" # "CAR" "independent"
source(paste(fldr,"scripts/philly_LOO.R", sep=""))
newseed_bool <- FALSE

############################### LOO & noSMOTE (Only CAR now) ###############################
smote_bool <- FALSE

############### PRETERM ###############
outcome <- "PRETERM" # "PRETERM","STILLBIRTH","CSECTION"

noRE_bool <- FALSE
prior <- "CAR" # "CAR" "independent"
source(paste(fldr,"scripts/philly_LOO.R", sep=""))



############### STILLBIRTH ############### 
outcome <- "STILLBIRTH" # "PRETERM","STILLBIRTH","CSECTION"

noRE_bool <- FALSE
prior <- "CAR" # "CAR" "independent"
source(paste(fldr,"scripts/philly_LOO.R", sep=""))
