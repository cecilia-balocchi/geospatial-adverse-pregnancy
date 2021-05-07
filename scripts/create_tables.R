library(xtable)
wdstr <- "results/"
##########################################################
############ Model comparison table (table 2) ############
##########################################################

add_str <- "noNH_WHITE_"; add_str2 <- "";

####### ALLDATA & noSMOTE
### PRETERM
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_",add_str,"newcov_PRETERM_CAR",add_str2,".RData"))
noSMOTE_CAR <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec), mean(WAIC1), mean(WAIC1_shifted), mean(WAIC2), mean(WAIC2_shifted))
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_",add_str,"newcov_PRETERM_independent",add_str2,".RData"))
noSMOTE_ind <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec), mean(WAIC1), mean(WAIC1_shifted), mean(WAIC2), mean(WAIC2_shifted))
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_",add_str,"noRE_newcov_PRETERM_independent",add_str2,".RData"))
noSMOTE_noRE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec), mean(WAIC1), mean(WAIC1_shifted), mean(WAIC2), mean(WAIC2_shifted))


res_pret_alld_noSMOTE <- rbind(noSMOTE_CAR,noSMOTE_ind, noSMOTE_noRE)
colnames(res_pret_alld_noSMOTE) <- c("AUC", "DIC","rmse","DIC*","rmse*","misclass", "sens", "spec","WAIC1","WAIC1*","WAIC2","WAIC2*")
res_pret_alld_noSMOTE <- res_pret_alld_noSMOTE[,c("DIC", "WAIC1", "WAIC2")]

### STILLBIRTH
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_",add_str,"newcov_STILLBIRTH_CAR",add_str2,".RData"))
noSMOTE_CAR <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec), mean(WAIC1), mean(WAIC1_shifted), mean(WAIC2), mean(WAIC2_shifted))
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_",add_str,"newcov_STILLBIRTH_independent",add_str2,".RData"))
noSMOTE_ind <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec), mean(WAIC1), mean(WAIC1_shifted), mean(WAIC2), mean(WAIC2_shifted))
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_",add_str,"noRE_newcov_STILLBIRTH_independent",add_str2,".RData"))
noSMOTE_noRE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec), mean(WAIC1), mean(WAIC1_shifted), mean(WAIC2), mean(WAIC2_shifted))


res_still_alld_noSMOTE <- rbind(noSMOTE_CAR,noSMOTE_ind, noSMOTE_noRE)
colnames(res_still_alld_noSMOTE) <- c("AUC", "DIC","rmse","DIC*","rmse*","misclass", "sens", "spec","WAIC1","WAIC1*","WAIC2","WAIC2*")
res_still_alld_noSMOTE <- res_still_alld_noSMOTE[,c("DIC", "WAIC1", "WAIC2")]


res <- t(rbind(res_still_alld_noSMOTE,res_pret_alld_noSMOTE))

xtable(res, digits = 1)

#############################################################################
############ Predictive performance, SMOTE vs ORIGINAL (table 3) ############
#############################################################################

method <- "independent"
####### LOO - PRETERM
tmp <-load(paste0(wdstr, "output_LOO_nogamma_rho_nointeractions_noNH_WHITE_newcov_PRETERM_",method,".RData"))
noSMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))
tmp <- load(paste0(wdstr,"output_LOO_nogamma_rho_SMOTE_nointeractions_noNH_WHITE_newcov_PRETERM_",method,".RData"))
SMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))

res_pret_loo <- rbind(SMOTE,noSMOTE)
colnames(res_pret_loo) <- c("AUC", "DIC","rmse","DIC*","rmse*","misclass", "sens", "spec")
res_pret_loo <- res_pret_loo[,-(2:5)]


####### ALLDATA - PRETERM
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_noNH_WHITE_newcov_PRETERM_",method,".RData"))
noSMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_SMOTE_nointeractions_noNH_WHITE_newcov_PRETERM_",method,".RData"))
SMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))

res_pret_all <- rbind(SMOTE,noSMOTE)
colnames(res_pret_all) <- c("AUC", "DIC","rmse","DIC*","rmse*","misclass", "sens", "spec")
res_pret_all <- res_pret_all[,-(2:5)]

res_pret <- rbind(res_pret_loo, res_pret_all)

####### LOO - STILLBIRTH
tmp <-load(paste0(wdstr,"output_LOO_nogamma_rho_nointeractions_noNH_WHITE_newcov_STILLBIRTH_CAR.RData"))
noSMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))
tmp <- load(paste0(wdstr,"output_LOO_nogamma_rho_SMOTE_nointeractions_noNH_WHITE_newcov_STILLBIRTH_CAR.RData"))
SMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))

res_still_loo <- rbind(SMOTE,noSMOTE)
colnames(res_still_loo) <- c("AUC", "DIC","rmse","DIC*","rmse*","misclass", "sens", "spec")
res_still_loo <- res_still_loo[,-(2:5)]


####### ALLDATA - STILLBIRTH
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_noNH_WHITE_newcov_STILLBIRTH_CAR.RData"))
noSMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_SMOTE_nointeractions_noNH_WHITE_newcov_STILLBIRTH_CAR.RData"))
SMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))

res_still_all <- rbind(SMOTE,noSMOTE)
colnames(res_still_all) <- c("AUC", "DIC","rmse","DIC*","rmse*","misclass", "sens", "spec")
res_still_all <- res_still_all[,-(2:5)]

res_still <- rbind(res_still_loo, res_still_all)

# round(t(res_pret),3)
# round(t(res_still),3)

res <- t(rbind(res_still, res_pret))
xtable(res, digits = 3)

#############################################################################
######################## Odds ratio table (table 4) #########################
#############################################################################

noNH_WHITE_bool <- TRUE
bayesp <- function(x){max(mean(x>0),mean(x<0))}
var_neighborhood <- c("prop_Asian", "prop_Hispanic", "prop_Black", # "Prop_White", 
                      "prop_women_15_to_50", "prop_women_below_poverty", 
                      "prop_women_public_assistance", "prop_women_labor_force", 
                      "prop_birth_last_12_months", "prop_women_HS_grad", 
                      "prop_women_college_grad", "log_occupied_housing", "log_housing_violation", 
                      "log_violent_crime", "log_nonviolent_crime")
if(noNH_WHITE_bool){
  var_individual <- c("Hispanic","Black", "Asian", "multiple_birth", "age")
} else {
  var_individual <- c("Hispanic","White","Black", "Asian", "multiple_birth", "age")
}


p1 <- length(var_individual)
p2 <- length(var_neighborhood)
col_to_be_scaled <- p1:(p1+p2)

burnin <- 500
thin <- 10
mcmc_niter <- 500*thin+burnin 
index_thinning <- seq(burnin, mcmc_niter, by = thin)

order_covariates <- c("age","Black","Hispanic","Asian","multiple_birth",
                      "prop_Asian","prop_Hispanic","prop_Black","prop_women_15_to_50",
                      "prop_women_below_poverty","prop_women_public_assistance",
                      "prop_women_labor_force","prop_birth_last_12_months",
                      "prop_women_HS_grad","prop_women_college_grad",
                      "log_occupied_housing","log_housing_violation",
                      "log_violent_crime","log_nonviolent_crime")
index <- match(order_covariates,c(var_individual, var_neighborhood))

### For PRETERM. Uncomment the lines with the stillbirth code to make the plot for STILLBIRTH.
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_noNH_WHITE_newcov_PRETERM_independent.RData"))
# tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_noNH_WHITE_newcov_STILLBIRTH_CAR.RData"))
year <- 1; output1 <- output_list1[[year]]; output2 <- output_list2[[year]]
beta_tr <- cbind(output1$beta[,index_thinning, drop = FALSE],
                 output2$beta[,index_thinning, drop = FALSE])
scale_sds <- scale_sds_list[[1]]
for(i in col_to_be_scaled){
  beta_tr[i,] <- beta_tr[i,] / scale_sds[i]
}
XtX <- XtX_list[[year]]
XtXinv <- solve(XtX)
VIF <- diag(XtX) * diag(XtXinv)
beta_OR <- rbind(exp(rowMeans(beta_tr[,index])),
                 exp(apply(beta_tr[,index], MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975))),
                 apply(beta_tr[,index], MARGIN = 1, FUN = bayesp),
                 VIF
)
colnames(beta_OR) <- c(var_individual, var_neighborhood)
rownames(beta_OR) <- c("Odds Ratio", "CI_LB", "CI_UB", "bayes_p", "VIF")
coeff_OR <- t(beta_OR) 
noSMOTE_OR <- coeff_OR[,c("Odds Ratio","bayes_p", "VIF")]

tmp <- load(paste0(wdstr,"output_alldata_nogamma_rho_SMOTE_nointeractions_noNH_WHITE_newcov_PRETERM_independent.RData"))
# tmp <- load(paste0(wdstr,"output_alldata_nogamma_rho_SMOTE_nointeractions_noNH_WHITE_newcov_STILLBIRTH_CAR.RData"))
year <- 1; output1 <- output_list1[[year]]; output2 <- output_list2[[year]]
beta_tr <- cbind(output1$beta[,index_thinning, drop = FALSE],
                 output2$beta[,index_thinning, drop = FALSE])
scale_sds <- scale_sds_list[[1]]
for(i in col_to_be_scaled){
  beta_tr[i,] <- beta_tr[i,] / scale_sds[i]
}
XtX <- XtX_list[[year]]
XtXinv <- solve(XtX)
VIF <- diag(XtX) * diag(XtXinv)
beta_OR <- rbind(exp(rowMeans(beta_tr[,index])),
                 exp(apply(beta_tr[,index], MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975))),
                 apply(beta_tr[,index], MARGIN = 1, FUN = bayesp),
                 VIF
)
colnames(beta_OR) <- c(var_individual, var_neighborhood)
rownames(beta_OR) <- c("Odds Ratio", "CI_LB", "CI_UB", "bayes_p", "VIF")
coeff_OR <- t(beta_OR) 
SMOTE_OR <- coeff_OR[,c("Odds Ratio","bayes_p", "VIF")]

rnames_noSM <- rownames(noSMOTE_OR)
rnames_SM <- rownames(SMOTE_OR)
res <- cbind(noSMOTE_OR, SMOTE_OR[match(rnames_noSM, rnames_SM),])
res_noVIF <- res[,-c(3,6)]
xtable(res_noVIF[index,])


#############################################################################
##################### Neighborhood comparison (in appendix) #################
#############################################################################

logit <- function(x) log(x/(1-x))
var_neighborhood <- c("prop_Asian", "prop_Hispanic", "prop_Black", # "Prop_White", 
                      "prop_women_15_to_50", "prop_women_below_poverty", 
                      "prop_women_public_assistance", "prop_women_labor_force", 
                      "prop_birth_last_12_months", "prop_women_HS_grad", 
                      "prop_women_college_grad", "log_occupied_housing", "log_housing_violation", 
                      "log_violent_crime", "log_nonviolent_crime")
var_neighborhood_old <- c("Prop_Asian", "Prop_Hispanic_Latino", "Prop_Black", #"Prop_White", 
                      "Prop_women15to50", "Prop_Below100percPoverty_women15to50", 
                      "Prop_ReceivedPublicAssistanceIncome_women15to50", "Prop_InLaborForce_women16to50", 
                      "Prop_BirthsInPast12Months_women15to50", "Prop_HighschoolGrad_women15to50", 
                      "Prop_BachelorsDegree_women15to50", "log.occupied.housing", "log.violation", 
                      "log.violent", "log.nonviolent")
if(noNH_WHITE_bool){
  var_individual <- c("Hispanic","Black", "Asian", "multiple_birth", "age")
} else {
  var_individual <- c("Hispanic","White","Black", "Asian", "multiple_birth", "age")
}
p1 <- length(var_individual)
p2 <- length(var_neighborhood)
col_to_be_scaled <- p1:(p1+p2)

burnin <- 500
thin <- 10
mcmc_niter <- 500*thin+burnin 
index_thinning <- seq(burnin, mcmc_niter, by = thin)

order_covariates <- c("age","Black","Hispanic","Asian","multiple_birth",
                      "prop_Asian","prop_Hispanic","prop_Black","prop_women_15_to_50",
                      "prop_women_below_poverty","prop_women_public_assistance",
                      "prop_women_labor_force","prop_birth_last_12_months",
                      "prop_women_HS_grad","prop_women_college_grad",
                      "log_occupied_housing","log_housing_violation",
                      "log_violent_crime","log_nonviolent_crime")
index <- match(order_covariates,c(var_individual, var_neighborhood))

philly_all_covariates$log.occupied.housing <- asinh(0.5*philly_all_covariates$Total_Occupied_Housing_Units)
philly_all_covariates$log.violation <- asinh(0.5*philly_all_covariates$violation)
philly_all_covariates$log.violent <- asinh(0.5*philly_all_covariates$violent)
philly_all_covariates$log.nonviolent <- asinh(0.5*philly_all_covariates$nonviolent)
# select only the covariates we use
philly_all_covariates <- (philly_all_covariates[,c(1,2,match(var_neighborhood_old, colnames(philly_all_covariates)))])

### we need this only for tractID_unique (tracts used in the analysis)
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_noNH_WHITE_newcov_STILLBIRTH_CAR.RData")) # STILL
tractID_unique <- tractID_unique_list[[1]] # used neighborhoods 
# select only the neighborhoods we use
philly_all_covariates <- philly_all_covariates[which(philly_all_covariates$tractID %in% tractID_unique),]


### find neighborhoods
IS1 <- c(); IS2 <- c(); IBS1 <- c(); IBS2 <- c()

####### Prop_black #########
# indeces for philly_all_covariates
i1 <- which.min(philly_all_covariates$Prop_Black); IS1 <- c(IS1,i1)
i2 <- which.max(philly_all_covariates$Prop_Black); IS2 <- c(IS2,i2)
# indeces for tractID_unique 
IBS1 <- c(IBS1,match(philly_all_covariates$tractID[i1],tractID_unique))
IBS2 <- c(IBS2,match(philly_all_covariates$tractID[i2],tractID_unique))
####### Prop_InLaborForce_women16to50 #######
# indeces for philly_all_covariates
i1 <- which.min(philly_all_covariates$Prop_InLaborForce_women16to50); IS1 <- c(IS1,i1)
i2 <- which.max(philly_all_covariates$Prop_InLaborForce_women16to50); IS2 <- c(IS2,i2)
# indeces for tractID_unique 
IBS1 <- c(IBS1,match(philly_all_covariates$tractID[i1],tractID_unique))
IBS2 <- c(IBS2,match(philly_all_covariates$tractID[i2],tractID_unique))
##########  Prop_BachelorsDegree_women15to50 #########
# indeces for philly_all_covariates
i1 <- which.min(philly_all_covariates$Prop_BachelorsDegree_women15to50); IS1 <- c(IS1,i1)
i2 <- which.max(philly_all_covariates$Prop_BachelorsDegree_women15to50); IS2 <- c(IS2,i2)
# indeces for tractID_unique 
IBS1 <- c(IBS1,match(philly_all_covariates$tractID[i1],tractID_unique))
IBS2 <- c(IBS2,match(philly_all_covariates$tractID[i2],tractID_unique))
###########  log.violent ###########
# indeces for philly_all_covariates
i1 <- which.min(philly_all_covariates$log.violent); IS1 <- c(IS1,i1)
i2 <- which.max(philly_all_covariates$log.violent); IS2 <- c(IS2,i2)
# indeces for tractID_unique 
IBS1 <- c(IBS1,match(philly_all_covariates$tractID[i1],tractID_unique))
IBS2 <- c(IBS2,match(philly_all_covariates$tractID[i2],tractID_unique))

# noSMOTE
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_noNH_WHITE_newcov_STILLBIRTH_CAR.RData")) # STILL
phat_neig <- phat_neig_list[[1]]

data <- data.frame()

i <- 1
i1 <- IS1[i]; i2 <- IS2[i]
iB1 <- IBS1[i]; iB2 <- IBS2[i]

data <- rbind(data, philly_all_covariates[c(i1,i2),])

data$phat[1:2] <- phat_neig[c(iB1,iB2)]
data$OR[1:2] <- exp(logit(phat_neig[c(iB1,iB2)]) )
data$RelOR[1:2] <- c(NA,exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(1,-1) )))

for(i in 2:4){
  i1 <- IS1[i]; i2 <- IS2[i]
  iB1 <- IBS1[i]; iB2 <- IBS2[i]
  data <- rbind(data, cbind(philly_all_covariates[c(i1,i2),], 
                            phat = phat_neig[c(iB1,iB2)],
                            OR = exp(logit(phat_neig[c(iB1,iB2)]) ),
                            RelOR = c(NA, exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(1,-1) )))))
}

data$phat2 <- NA
data$OR2 <- NA
data$RelOR2 <- NA

tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_noNH_WHITE_newcov_PRETERM_independent.RData")) # PRE
phat_neig <- phat_neig_list[[1]]
for(i in 1:4){
  i1 <- IS1[i]; i2 <- IS2[i]
  iB1 <- IBS1[i]; iB2 <- IBS2[i]
  
  data$phat2[2*(i-1)+1:2] = phat_neig[c(iB1,iB2)]
  data$OR2[2*(i-1)+1:2] = exp(logit(phat_neig[c(iB1,iB2)]) )
  data$RelOR2[2*(i-1)+1:2] = c(NA, exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(1,-1) )))
}


data_print <- data
data_print$tractID <- as.character(data_print$tractID)
colnames(data_print)[match(var_neighborhood_old,colnames(data_print))] <- var_neighborhood
data_print[, var_neighborhood] <- round(data_print[, var_neighborhood],2)
data_print[, c("phat","OR","RelOR","phat2","OR2","RelOR2")] <- round(data_print[, c("phat","OR","RelOR","phat2","OR2","RelOR2")],3)

xtable(t(data_print))
