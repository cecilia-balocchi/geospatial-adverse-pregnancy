binary_to_factor <- function(dataframe){
  ## dataframe needs to have columns c("NH_AFAM","NH_WHITE","NH_ASIAN","HISPANIC")
  ## returns a dataframe with the same columns + column RACE
  dataframe$RACE <- NULL
  dataframe$RACE[dataframe[,"NH_AFAM"] == 1] <- "NH_AFAM"
  dataframe$RACE[dataframe[,"NH_WHITE"] == 1] <- "NH_WHITE"
  dataframe$RACE[dataframe[,"NH_ASIAN"] == 1] <- "NH_ASIAN"
  dataframe$RACE[dataframe[,"HISPANIC"] == 1] <- "HISPANIC"
  dataframe$RACE <- as.factor(dataframe$RACE)
  return(dataframe)
}
factor_to_binary <- function(dataframe){
  ## dataframe needs to have column RACE 
  ## it can also have c("NH_AFAM","NH_WHITE","NH_ASIAN","HISPANIC"), or not
  dataframe$NH_AFAM <- 0
  dataframe$NH_WHITE <- 0
  dataframe$NH_ASIAN <- 0
  dataframe$HISPANIC <- 0
  dataframe$NH_AFAM[dataframe$RACE == "NH_AFAM"] <- 1
  dataframe$NH_WHITE[dataframe$RACE == "NH_WHITE"] <- 1
  dataframe$NH_ASIAN[dataframe$RACE == "NH_ASIAN"] <- 1
  dataframe$HISPANIC[dataframe$RACE == "HISPANIC"] <- 1
  return(dataframe)
}