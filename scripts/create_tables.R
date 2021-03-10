##########################################################
############ Model comparison table (table 2) ############
##########################################################

library(xtable)
setwd("../cecis_scripts/")

wdstr <- "results/results_feb8_9/"
####### ALLDATA 
### noSMOTE
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_newcov_PRETERM_CAR_correct.RData"))
noSMOTE_CAR <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec), mean(WAIC1), mean(WAIC1_shifted), mean(WAIC2), mean(WAIC2_shifted))
rho_noSMOTE <- quantile(output_list[[1]]$rho, probs = c(0.025,0.25,0.5,0.75,0.975))
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_newcov_PRETERM_independent_correct.RData"))
noSMOTE_ind <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec), mean(WAIC1), mean(WAIC1_shifted), mean(WAIC2), mean(WAIC2_shifted))
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_noRE_newcov_PRETERM_independent_correct.RData"))
noSMOTE_noRE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec), mean(WAIC1), mean(WAIC1_shifted), mean(WAIC2), mean(WAIC2_shifted))


res_pret_alld_noSMOTE <- rbind(noSMOTE_CAR,noSMOTE_ind, noSMOTE_noRE)
colnames(res_pret_alld_noSMOTE) <- c("AUC", "DIC","rmse","DIC*","rmse*","misclass", "sens", "spec","WAIC1","WAIC1*","WAIC2","WAIC2*")
res_pret_alld_noSMOTE <- res_pret_alld_noSMOTE[,c("DIC", "WAIC1", "WAIC2")]

### noSMOTE
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_newcov_STILLBIRTH_CAR_correct.RData"))
noSMOTE_CAR <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec), mean(WAIC1), mean(WAIC1_shifted), mean(WAIC2), mean(WAIC2_shifted))
rho_noSMOTE <- quantile(output_list[[1]]$rho, probs = c(0.025,0.25,0.5,0.75,0.975))
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_newcov_STILLBIRTH_independent_correct.RData"))
noSMOTE_ind <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec), mean(WAIC1), mean(WAIC1_shifted), mean(WAIC2), mean(WAIC2_shifted))
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_noRE_newcov_STILLBIRTH_independent_correct.RData"))
noSMOTE_noRE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec), mean(WAIC1), mean(WAIC1_shifted), mean(WAIC2), mean(WAIC2_shifted))


res_still_alld_noSMOTE <- rbind(noSMOTE_CAR,noSMOTE_ind, noSMOTE_noRE)
colnames(res_still_alld_noSMOTE) <- c("AUC", "DIC","rmse","DIC*","rmse*","misclass", "sens", "spec","WAIC1","WAIC1*","WAIC2","WAIC2*")
res_still_alld_noSMOTE <- res_still_alld_noSMOTE[,c("DIC", "WAIC1", "WAIC2")]


res <- t(rbind(res_still_alld_noSMOTE,res_pret_alld_noSMOTE))

xtable(res, digits = 1)

#############################################################################
############ Predictive performance, SMOTE vs ORIGINAL (table 3) ############
#############################################################################

library(xtable)
setwd("../cecis_scripts/")

wdstr <- "results/results_feb8_9/"

####### LOO - PRETERM
tmp <-load("results/Cecilia_code_output_Dec29/output_LOO_nogamma_rho_nointeractions_newcov_PRETERM_CAR.RData")
noSMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))
tmp <- load("results/Cecilia_code_output_Dec29/output_LOO_nogamma_rho_SMOTE_nointeractions_newcov_PRETERM_CAR.RData")
SMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))

res_pret_loo <- rbind(SMOTE,noSMOTE)
colnames(res_pret_loo) <- c("AUC", "DIC","rmse","DIC*","rmse*","misclass", "sens", "spec")
res_pret_loo <- res_pret_loo[,-(2:5)]


####### ALLDATA - PRETERM
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_newcov_PRETERM_CAR_correct.RData"))
noSMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_SMOTE_nointeractions_newcov_PRETERM_CAR_correct.RData"))
SMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))

res_pret_all <- rbind(SMOTE,noSMOTE)
colnames(res_pret_all) <- c("AUC", "DIC","rmse","DIC*","rmse*","misclass", "sens", "spec")
res_pret_all <- res_pret_all[,-(2:5)]

res_pret <- rbind(res_pret_loo, res_pret_all)

####### LOO - STILLBIRTH
tmp <-load("results/Cecilia_code_output_Dec29/output_LOO_nogamma_rho_nointeractions_newcov_STILLBIRTH_CAR.RData")
noSMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))
tmp <- load("results/Cecilia_code_output_Dec29/output_LOO_nogamma_rho_SMOTE_nointeractions_newcov_STILLBIRTH_CAR.RData")
SMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))

res_still_loo <- rbind(SMOTE,noSMOTE)
colnames(res_still_loo) <- c("AUC", "DIC","rmse","DIC*","rmse*","misclass", "sens", "spec")
res_still_loo <- res_still_loo[,-(2:5)]


####### ALLDATA - STILLBIRTH
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_newcov_STILLBIRTH_CAR_correct.RData"))
noSMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_SMOTE_nointeractions_newcov_STILLBIRTH_CAR_correct.RData"))
SMOTE <- c(mean(auc), mean(DIC), mean(rmse), mean(DIC_shifted), mean(rmse_shifted), mean(misclass), mean(sens), mean(spec))

res_still_all <- rbind(SMOTE,noSMOTE)
colnames(res_still_all) <- c("AUC", "DIC","rmse","DIC*","rmse*","misclass", "sens", "spec")
res_still_all <- res_still_all[,-(2:5)]

res_still <- rbind(res_still_loo, res_still_all)

round(t(res_pret),4)
round(t(res_still),4)

res <- t(rbind(res_still, res_pret))
xtable(res, digits = 4)

############ Odds ratio tables (table 4 and 5) ############

setwd("../cecis_scripts/")

noMULTIPLE_bool <- FALSE
bayesp <- function(x){max(mean(x>0),mean(x<0))}
var_neighborhood <- c("Prop_Asian", "Prop_Hispanic", "Prop_Black", # "Prop_White", 
                      "prop_women_15_to_50", "prop_women_below_poverty", 
                      "prop_women_public_assistance", "prop_women_labor_force", 
                      "prop_birth_last_12_months", "prop_women_HS_grad", 
                      "prop_women_college_grad", "log_occupied_housing", "log_housing_violation", 
                      "log_violent_crime", "log_nonviolent_crime")
if(noMULTIPLE_bool){
  var_individual <- c("Hispanic","white", "black", "Asian", "age")
} else {
  var_individual <- c("Hispanic","white", "black", "Asian", "multiple_birth", 
                      "age")
}

p1 <- length(var_individual)
p2 <- length(var_neighborhood)
col_to_be_scaled <- p1:(p1+p2)

burnin <- 500
thin <- 5
mcmc_niter <- 2001
index <- seq(burnin, mcmc_niter, by = thin)

order_covariates <- c("age","white","black","Hispanic","Asian","multiple_birth",
                      "Prop_Asian","Prop_Hispanic","Prop_Black","prop_women_15_to_50",
                      "prop_women_below_poverty","prop_women_public_assistance",
                      "prop_women_labor_force","prop_birth_last_12_months",
                      "prop_women_HS_grad","prop_women_college_grad",
                      "log_occupied_housing","log_housing_violation",
                      "log_violent_crime","log_nonviolent_crime")
index <- match(order_covariates,c(var_individual, var_neighborhood))
### PRE
# no white
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_newcov_PRETERM_CAR_correct.RData"))
tmp <- load(paste0(wdstr,"output_alldata_nogamma_rho_SMOTE_nointeractions_newcov_PRETERM_CAR_correct.RData"))

### STILL
# no white
tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_newcov_STILLBIRTH_CAR_correct.RData"))
tmp <- load(paste0(wdstr,"output_alldata_nogamma_rho_SMOTE_nointeractions_newcov_STILLBIRTH_CAR_correct.RData"))


year <- 1
output <- output_list[[year]]

beta_tr <- output$beta
scale_sds <- scale_sds_list[[1]]
for(i in col_to_be_scaled){
  beta_tr[i,] <- beta_tr[i,] / scale_sds[i]
}

XtX <- XtX_list[[year]]
# n <- sum(diag(XtX[1:4,1:4]))-1
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
# coeff_OR <- coeff_OR[sort(coeff_OR[,"Odds Ratio"], index.return = T, decreasing = T)$ix,]


noSMOTE_OR <- coeff_OR[,c("Odds Ratio","bayes_p", "VIF")]
SMOTE_OR <- coeff_OR[,c("Odds Ratio","bayes_p", "VIF")]

rnames_noSM <- rownames(noSMOTE_OR)
rnames_SM <- rownames(SMOTE_OR)

res <- cbind(noSMOTE_OR, SMOTE_OR[match(rnames_noSM, rnames_SM),])
xtable(res)

res_noVIF <- res[,-c(3,6)]
xtable(res_noVIF[index,])