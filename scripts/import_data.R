setwd("~/Google Drive/My Drive/STAT/Research/3Urban/Maternal_health/code/")
## read fake datasets
setwd("Geographic_analysis/")
# encounter_at_deliv_woutcomes_trimmed <- read.csv("RANDOM_DATA_AND_SCRIPT/encounter_at_deliv_woutcomes_trimmed_RANDOM.csv")

# note: rownames are different but who cares
data.new.philly <- read.csv("RANDOM_DATA_AND_SCRIPT/data.new.philly_RANDOM_NEW.csv")
# data.new.philly$tractID <- as.character(data.new.philly$tractID)

# counts_outcomes_philly_geo <- read.csv("RANDOM_DATA_AND_SCRIPT/counts_outcomes_philly_geo_RANDOM.csv")
# counts_outcomes_philly_geo <- counts_outcomes_philly_geo[,-1]
# counts_outcomes_philly_geo$tractID <- as.factor(as.character(counts_outcomes_philly_geo$tractID))

## read other datasets that were imported in script1

data_path <- "data/"
load(paste0(data_path,"phillytracts"))



library(spdep)

tract_df = tracts@data

## Fake data. Uncomment this out if you are working on the fake data
## data.new.philly <- read.csv("data.new.philly_RANDOM.csv",header=TRUE,sep=",")

## Load covariates data
# fldr="/Users/bolandm/Box\ Sync/Geographic_analysis/"
fldr=""

philly_all_covariates <- read.csv(paste(fldr, "data/ALL_NEIGHBORHOOD_COVARIATES_YEARLY_PHILLY.csv",sep=""),
                                        header=TRUE,sep=",")
# head(philly_all_covariates)
# colnames(philly_all_covariates)
## Remove 1st, 44th, and 45 columns
#philly_all_covariates <- philly_all_covariates[,-c(1,45:46)]  #i changed this on 8/7/2020 removed imm_dang and unsafe
philly_all_covariates <- philly_all_covariates[,-c(48:49)] ## Remove imm_dang and unsafe as too many NA's

## Change column name "GEOID10" to "tractID"
colnames(philly_all_covariates)[1] <- c("tractID")
# colnames(philly_all_covariates)


### careful!! tract is NOT tractID
if(!all(data.new.philly[,14] < 400)){
  colnames(data.new.philly)[14]<- c("tractID")
}
# colnames(data.new.philly)[14]<- c("tractID")
# philly_all_covariates$tractID <- as.character(philly_all_covariates$tractID)

## Merge covariates with data.new.philly
data.new.philly <- base::merge(data.new.philly, philly_all_covariates, by=c("tractID","YEAR"), all=FALSE)
data.new.philly$tractID <- as.character(data.new.philly$tractID)

## Number of total observations, total number of preterm births and stillbirths 
total_tracts = length(unique(data.new.philly$tractID))
total_observations = dim(data.new.philly)[1]
preterm_col = data.new.philly[,colnames(data.new.philly)=="PRETERM"]
total_preterm = length(which(preterm_col==1))
stillbirth_col = data.new.philly[,colnames(data.new.philly)=="STILLBIRTH"]
total_stillbirth = length(which(stillbirth_col==1))

## Let's select which tracts to keep, considering the leave-one-out and that we want to have
## at least 10 deliveries in each tract in each iteration of the leave-one-out 
# tmp = number of deliveries in each tract for each year
# tmp2 = is the number of deliveries that we can use with loo when we exclude that year
# [tmp2 is initialized to the array with total number of deliveries 
#  in each tract (over all the years) (same number per row),
#  then we subtract the number of deliveries in that year]
# tmp3 = TRUE if in one year there are less than 10 deliveries
tmp = table(data.new.philly$tractID, data.new.philly$YEAR)
tmp2 = array(data = rep(rowSums(tmp), 8), dim = c(384,8))
for(i in 1:384){
  tmp2[i,] = tmp2[i,] - tmp[i,]
}
tmp3 = (tmp2<10)

## ind_atleast10c are the tracts with more than 10 deliveries in all loo years (NO loo years with less than 10)
ind_atleast10c = which(apply(tmp3, MARGIN = 1, function(x)ifelse(any(x), 1, 0))==0)
# length(ind_atleast10c) # 296 tracts - only 88 removed 
tractID_tokeep = rownames(tmp)[ind_atleast10c]

# let's check that the tracts with NAs in the covariates are excluded:
var_neighborhood1 <- c("Prop_Asian", "Prop_Hispanic_Latino", "Prop_Black", 
                      "Prop_White", "Prop_women15to50", "Prop_Below100percPoverty_women15to50", 
                      "Prop_ReceivedPublicAssistanceIncome_women15to50", "Prop_InLaborForce_women16to50", 
                      "Prop_BirthsInPast12Months_women15to50", "Prop_HighschoolGrad_women15to50", 
                      "Prop_BachelorsDegree_women15to50", "Total_Occupied_Housing_Units", "violation", 
                      "violent", "nonviolent")
index.data.NA <- which(apply(X = data.new.philly[,var_neighborhood1], MARGIN = 1, FUN = function(x) any(is.na(x))))
tractID.NA <- unique(data.new.philly[index.data.NA,"tractID"])
# tractID.after.NA <- unique(data.new.philly[-index.data.NA,"tractID"])
# tractID_tokeep <- tractID_tokeep[which(tractID_tokeep %in% tractID.after.NA)]

## remember that the new covariates vary over time, so some years can be NA and others no.
## so, subsetting based on the neighbhoords will leave some NAs!
## problem is: it's possible when we exclude 1 year, we get less neighborhoods that the ones here
## so now instead we exclude the neighborhoods that have NA in at least one year.

tractID_tokeep <- setdiff(tractID_tokeep, tractID.NA)

tractID_tokeep1 = tractID_tokeep # to get restricted w.sym


############################################################################
############################################################################
## 
## IF WE WANT TO RUN CAR WE NEED TO RUN THIS PART 
## ALL NEIGHBORHOODS SHOULD BE CONNECTED FOR CAR
## SO LET'S REMOVE ISOLATED NEIGHBORHOODS

# load(paste(fldr,"data/phillytracts", sep=""))
# tract_df = tracts@data

temp = tract_df$GEOID10 %in% tractID_tokeep1
index.wsym = as.numeric(rownames(tract_df[temp,]))+1
# index.wsym = (1:384)[temp] # this is equivalent

n_tr <- 384
list.poly <- poly2nb(tracts)
w <- matrix(data = 0, nrow = n_tr, ncol = n_tr)
for(i in 1:n_tr){
  w[i, list.poly[[i]]] <- 1
}
w.sym <- w + t(w) - w*t(w)
w.sym <- w.sym[index.wsym, index.wsym]
which(rowSums(w.sym)==0)
# if something is returned, it means some rows are indeed isolated and we need to remove them

if( sum(rowSums(w.sym)==0) > 0){
  ## Since some are isolated, let's remove them
  exclude_index = index.wsym[which(rowSums(w.sym) == 0)] # these are the rows of tract_df that we should also exclude
  ## check that they were not already excluded
  levels(tract_df$GEOID10)[tract_df$GEOID10[exclude_index]] %in% tractID_tokeep # TRUE
  
  tractID_tokeep = tractID_tokeep1[-match(levels(tract_df$GEOID10)[tract_df$GEOID10[exclude_index]], tractID_tokeep)]
  temp = tract_df$GEOID10 %in% tractID_tokeep
  index.wsym = as.numeric(rownames(tract_df[temp,]))+1
  w.sym <- w + t(w) - w*t(w)
  w.sym <- w.sym[index.wsym, index.wsym]
  which(rowSums(w.sym)==0) # should now be empty
}

############################################################################
## LET'S SAVE THE NEIGHBORHOODS WE KEEP!
# choose the right path and then upload on box 
# (you can overwrite the ones that were there before)

# save(tractID_tokeep, file = paste0(fldr,"data/tractID_tokeep.Rdata"))
# save(index.wsym, file = paste0(fldr,"data/index_wsym.Rdata"))
############################################################################


dim(data.new.philly)
index_tokeep = which(data.new.philly$tractID %in% tractID_tokeep)
## this is data.new.philly restricted to the tracts with at least 10 - USE THIS
data.philly.atleast10 = data.new.philly[index_tokeep,]

## total observations, total number of preterm births and stillbirths 
## after preprocessing (i.e. removing missing data and tracts with <10 pregnancies or isolated tracts) 
tracts_after_preprocessing = length(tractID_tokeep)
total_observations_after_preprocessing = dim(data.philly.atleast10)[1]
preterm_after_preprocessing = data.philly.atleast10$PRETERM
total_preterm_after_preprocessing = length(which(preterm_after_preprocessing==1))
stillbirth_after_preprocessing = data.philly.atleast10$STILLBIRTH
total_stillbirth_after_preprocessing = length(which(stillbirth_after_preprocessing==1))

## Transform the number of occupied housing units,
## the umber of housing violations, 
## the number of violent crimes,
## and nonviolent crimes with a robust log transform
data.philly.atleast10$log.occupied.housing <- asinh(0.5*data.philly.atleast10$Total_Occupied_Housing_Units)
data.philly.atleast10$log.violation <- asinh(0.5*data.philly.atleast10$violation)
data.philly.atleast10$log.violent <- asinh(0.5*data.philly.atleast10$violent)
data.philly.atleast10$log.nonviolent <- asinh(0.5*data.philly.atleast10$nonviolent)

#making race variable numeric
data.philly.atleast10$NH_AFAM = as.numeric(as.character(data.philly.atleast10$NH_AFAM))
data.philly.atleast10$NH_WHITE = as.numeric(as.character(data.philly.atleast10$NH_WHITE))
data.philly.atleast10$NH_ASIAN = as.numeric(as.character(data.philly.atleast10$NH_ASIAN))
data.philly.atleast10$HISPANIC = as.numeric(as.character(data.philly.atleast10$HISPANIC))

## Include interaction terms for race
data.philly.atleast10$black_interaction <- data.philly.atleast10$NH_AFAM*data.philly.atleast10$Prop_Black
data.philly.atleast10$white_interaction <- data.philly.atleast10$NH_WHITE*data.philly.atleast10$Prop_White
data.philly.atleast10$asian_interaction <- data.philly.atleast10$NH_ASIAN*data.philly.atleast10$Prop_Asian
data.philly.atleast10$hispanic_interaction <- data.philly.atleast10$HISPANIC*data.philly.atleast10$Prop_Hispanic_Latino


# data.preg.tractcov.atleast10 columns that we care about
data.preg.tractcov.atleast10 = data.philly.atleast10[,c("tractID", "YEAR", "PREG_ID",
                                                        "LONGITUDE", "LATITUDE",
                                                        "NH_AFAM", "NH_WHITE", "NH_ASIAN", "HISPANIC",
                                                        "ENC_AGE_DEC", "MULTIPLE_BIRTH", "STILLBIRTH", "PRETERM", "CSECTION",
                                                        "Prop_Asian", "Prop_Hispanic_Latino", "Prop_Black", "Prop_White",
                                                        "Prop_women15to50", "Prop_Below100percPoverty_women15to50",
                                                        "Prop_ReceivedPublicAssistanceIncome_women15to50",
                                                        "Prop_InLaborForce_women16to50", "Prop_BirthsInPast12Months_women15to50",
                                                        "Prop_HighschoolGrad_women15to50", "Prop_BachelorsDegree_women15to50",
                                                        "log.occupied.housing", "log.violation", "log.violent", "log.nonviolent",
                                                        "black_interaction", "white_interaction", "asian_interaction", "hispanic_interaction")]



# head(data.preg.tractcov.atleast10)
# dim(data.preg.tractcov.atleast10)
#45919 pregnancies!

# PRETERM IS A FACTOR!! [NOT with the fake data!!]
# data.preg.tractcov.atleast10$PRETERM = as.numeric(levels(data.preg.tractcov.atleast10$PRETERM))[data.preg.tractcov.atleast10$PRETERM]

# 
# ## Export a CSV file
# before.after.numbers <- c(total_tracts,
#                           total_observations,
#                           total_preterm,
#                           total_stillbirth,
#                           tracts_after_preprocessing,
#                           total_observations_after_preprocessing,
#                           total_preterm_after_preprocessing,
#                           total_stillbirth_after_preprocessing)
# 
# names(before.after.numbers) <- c("Total Tracts Before Preprocessing",
#                                  "Total Observations Before Preprocessing",
#                                  "Total Preterm Before Preprocessing",
#                                  "Total Stillbirth Before Preprocessing",
#                                  "Total Tracts After Preprocessing",
#                                  "Total Observations After Preprocessing",
#                                  "Total Preterm After Preprocessing",
#                                  "Total Stillbirth After Preprocessing")
# before.after.numbers <- as.data.frame(before.after.numbers)

# Write CSV file
# write.csv(before.after.numbers, paste(fldr, "scripts/before_and_after_summary_statistics.csv",sep=""))

