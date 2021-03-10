setwd("../cecis_scripts/")
wdstr <- "results/results_feb8_9/"
logit <- function(x) log(x/(1-x))

noMULTIPLE_bool <- FALSE
bayesp <- function(x){max(mean(x>0),mean(x<0))}
var_neighborhood <- c("Prop_Asian", "Prop_Hispanic_Latino", "Prop_Black", #"Prop_White",
                      "Prop_women15to50", "Prop_Below100percPoverty_women15to50", 
                      "Prop_ReceivedPublicAssistanceIncome_women15to50", "Prop_InLaborForce_women16to50", 
                      "Prop_BirthsInPast12Months_women15to50", "Prop_HighschoolGrad_women15to50", 
                      "Prop_BachelorsDegree_women15to50", "log.occupied.housing", "log.violation", 
                      "log.violent", "log.nonviolent")
if(noMULTIPLE_bool){
  var_individual <- c("HISPANIC","NH_WHITE", "NH_AFAM", "NH_ASIAN", "ENC_AGE_DEC")
} else {
  var_individual <- c("HISPANIC","NH_WHITE", "NH_AFAM", "NH_ASIAN", "MULTIPLE_BIRTH", 
                      "ENC_AGE_DEC")
}

p1 <- length(var_individual)
p2 <- length(var_neighborhood)
col_to_be_scaled <- p1:(p1+p2)

burnin <- 500
thin <- 10 #5
mcmc_niter <- 2001
index <- seq(burnin, mcmc_niter, by = thin)

order_covariates <- c("age","white","black","Hispanic","Asian","multiple_birth",
                      "Prop_Asian","Prop_Hispanic","Prop_Black",
                      "prop_women_15_to_50","prop_women_below_poverty","prop_women_public_assistance",
                      "prop_women_labor_force","prop_birth_last_12_months",
                      "prop_women_HS_grad","prop_women_college_grad",
                      "log_occupied_housing","log_housing_violation",
                      "log_violent_crime","log_nonviolent_crime")
# index <- match(order_covariates,c(var_individual, var_neighborhood))

philly_all_covariates$log.occupied.housing <- asinh(0.5*philly_all_covariates$Total_Occupied_Housing_Units)
philly_all_covariates$log.violation <- asinh(0.5*philly_all_covariates$violation)
philly_all_covariates$log.violent <- asinh(0.5*philly_all_covariates$violent)
philly_all_covariates$log.nonviolent <- asinh(0.5*philly_all_covariates$nonviolent)
# select only the covariates we use
philly_all_covariates <- (philly_all_covariates[,c(1,2,match(var_neighborhood, colnames(philly_all_covariates)))])


# PRE, noSMOTE
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_newcov_PRETERM_CAR_correct.RData"))


year <- 1
output <- output_list[[year]]

alphas <- rowMeans(output$alpha[,index])
beta_tr <- output$beta
scale_sds <- scale_sds_list[[1]]
for(i in col_to_be_scaled){
  beta_tr[i,] <- beta_tr[i,] / scale_sds[i]
}
beta_LOR <- rbind(rowMeans(beta_tr[,index]),
                  apply(beta_tr[,index], MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975)))
colnames(beta_LOR) <- c(var_individual, var_neighborhood)
rownames(beta_LOR) <- c("Log Odds Ratio", "CI_LB", "CI_UB")
beta_tr_neigh <- beta_LOR[1, var_neighborhood]
# used neighborhoods 
tractID_unique <- tractID_unique_list[[1]]
# select only the neighborhoods we use
philly_all_covariates <- philly_all_covariates[which(philly_all_covariates$tractID %in% tractID_unique),]

IQRs <- apply(philly_all_covariates[,var_neighborhood], MARGIN = 2, IQR, na.rm = T)

phat_neig <- phat_neig_list[[1]]


data <- data.frame()

######################################## Prop_black ########################################
# indeces for philly_all_covariates
i1 <- which.min(philly_all_covariates$Prop_Black) # 81
i2 <- which.max(philly_all_covariates$Prop_Black) # 883
# indeces for tractID_unique 
iB1 <- match(philly_all_covariates$tractID[i1],tractID_unique)
iB2 <- match(philly_all_covariates$tractID[i2],tractID_unique)

t(philly_all_covariates[c(i1,i2),])

data <- rbind(data, philly_all_covariates[c(i1,i2),])
data$OR <- NA
data$RelOR <- NA

abs(philly_all_covariates[i1,var_neighborhood]-philly_all_covariates[i2,var_neighborhood])/IQRs

alphas[match(philly_all_covariates$tractID[i1],tractID_unique)]

alpha_diff <- alphas[iB1] - alphas[iB2] 
beta_diff <- as.numeric(philly_all_covariates[i1,var_neighborhood]-philly_all_covariates[i2,var_neighborhood])
exp(alpha_diff + beta_diff %*% beta_tr_neigh)
exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(1,-1) ))
data$phat[1:2] <- phat_neig[c(iB1,iB2)]
data$OR[1:2] <- exp(logit(phat_neig[c(iB1,iB2)]) )
data$RelOR[1:2] <- c(NA,exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(1,-1) )))



######################################## Prop_InLaborForce_women16to50 ######################################## 
# indeces for philly_all_covariates
i1 <- which.min(philly_all_covariates$Prop_InLaborForce_women16to50) 
i2 <- which.max(philly_all_covariates$Prop_InLaborForce_women16to50) 
# indeces for tractID_unique 
iB1 <- match(philly_all_covariates$tractID[i1],tractID_unique)
iB2 <- match(philly_all_covariates$tractID[i2],tractID_unique)

# t(philly_all_covariates[c(i1,i2),])
# abs(philly_all_covariates[i1,var_neighborhood]-philly_all_covariates[i2,var_neighborhood])/IQRs

alpha_diff <- alphas[iB1] - alphas[iB2] 
beta_diff <- as.numeric(philly_all_covariates[i1,var_neighborhood]-philly_all_covariates[i2,var_neighborhood])
exp(alpha_diff + beta_diff %*% beta_tr_neigh)
exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(1,-1) ))

data <- rbind(data, cbind(philly_all_covariates[c(i1,i2),], 
                          phat = phat_neig[c(iB1,iB2)],
                          OR = exp(logit(phat_neig[c(iB1,iB2)]) ),
                          RelOR = c(NA, exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(1,-1) )))))


########################################  Prop_BachelorsDegree_women15to50 ######################################## 
# indeces for philly_all_covariates
i1 <- which.min(philly_all_covariates$Prop_BachelorsDegree_women15to50) 
i2 <- which.max(philly_all_covariates$Prop_BachelorsDegree_women15to50) 
# indeces for tractID_unique 
iB1 <- match(philly_all_covariates$tractID[i1],tractID_unique)
iB2 <- match(philly_all_covariates$tractID[i2],tractID_unique)

# t(philly_all_covariates[c(i1,i2),])
# abs(philly_all_covariates[i1,var_neighborhood]-philly_all_covariates[i2,var_neighborhood])/IQRs

alpha_diff <- alphas[iB1] - alphas[iB2] 
beta_diff <- as.numeric(philly_all_covariates[i1,var_neighborhood]-philly_all_covariates[i2,var_neighborhood])
exp(alpha_diff + beta_diff %*% beta_tr_neigh)
exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(1,-1) ))

data <- rbind(data, cbind(philly_all_covariates[c(i1,i2),], 
                          phat = phat_neig[c(iB1,iB2)],
                          OR = exp(logit(phat_neig[c(iB1,iB2)]) ),
                          RelOR = c(NA, exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(1,-1) )))))

########################################  log.violent ######################################## 
# indeces for philly_all_covariates
i1 <- which.min(philly_all_covariates$log.violent) 
i2 <- which.max(philly_all_covariates$log.violent) 
# indeces for tractID_unique 
iB1 <- match(philly_all_covariates$tractID[i1],tractID_unique)
iB2 <- match(philly_all_covariates$tractID[i2],tractID_unique)

# t(philly_all_covariates[c(i1,i2),])
# abs(philly_all_covariates[i1,var_neighborhood]-philly_all_covariates[i2,var_neighborhood])/IQRs

alpha_diff <- alphas[iB1] - alphas[iB2] 
beta_diff <- as.numeric(philly_all_covariates[i1,var_neighborhood]-philly_all_covariates[i2,var_neighborhood])
exp(alpha_diff + beta_diff %*% beta_tr_neigh)
exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(1,-1) ))

data <- rbind(data, cbind(philly_all_covariates[c(i1,i2),], 
                          phat = phat_neig[c(iB1,iB2)],
                          OR = exp(logit(phat_neig[c(iB1,iB2)]) ),
                          RelOR = c(NA, exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(1,-1) )))))



data_print <- data
data_print$tractID <- as.character(data_print$tractID)
data_print[, var_neighborhood] <- round(data_print[, var_neighborhood],2)
data_print[, c("phat","OR","RelOR")] <- round(data_print[, c("phat","OR","RelOR")],3)

colnames(data_print)[3:(2+length(var_neighborhood))] <- order_covariates[7:length(order_covariates)]


library(xtable)
xtable(t(data_print))
