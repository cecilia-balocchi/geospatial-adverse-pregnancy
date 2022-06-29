library(xtable)
wdstr <- "results/"

#############################################################################
######################## Odds ratio table (table 2) #########################
#############################################################################

bayesp <- function(x){max(mean(x>0),mean(x<0))}
var_neighborhood <- c("proportion Asian","proportion Hispanic","proportion Black", 
                      "proportion women", "poverty", 
                      "public assistance", "labor force", 
                      "recent birth", "high school grad", 
                      "college grad", "occupied housing", "housing violation", 
                      "violent crime", "nonviolent crime")
var_individual <- c("Hispanic","Black", "Asian", "multiple birth", "age")

p1 <- length(var_individual)
p2 <- length(var_neighborhood)
col_to_be_scaled <- p1:(p1+p2)

burnin <- 500
thin <- 5#10
mcmc_niter <- 1000*thin+burnin
index_thinning <- seq(burnin, mcmc_niter, by = thin)

order_covariates <- c("age","Black","Hispanic","Asian","multiple birth",
                      "proportion Asian","proportion Hispanic","proportion Black","proportion women",
                      "poverty","public assistance",
                      "labor force","recent birth",
                      "high school grad","college grad",
                      "occupied housing","housing violation",
                      "violent crime","nonviolent crime")
index <- match(order_covariates,c(var_individual, var_neighborhood))

year <- 8; add_str = ""; add_str2 <- "_YEAR8"
for(output_string in c("PRETERM", "STILLBIRTH")){
  
  tmp <-load(paste0(wdstr,"output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
  output1 <- output_list1[[year]]; output2 <- output_list2[[year]]
  beta_tr <- cbind(output1$beta[,index_thinning, drop = FALSE],
                   output2$beta[,index_thinning, drop = FALSE])
  scale_sds <- scale_sds_list[[year]]
  for(i in col_to_be_scaled){
    beta_tr[i,] <- beta_tr[i,] / scale_sds[i]
  }
  XtX <- XtX_list[[year]]
  XtXinv <- solve(XtX)
  VIF <- diag(XtX) * diag(XtXinv)
  beta_OR <- rbind(exp(rowMeans(beta_tr[,index])),
                   exp(apply(beta_tr[,index], MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975))),
                   apply(beta_tr[,index], MARGIN = 1, FUN = bayesp),
                   VIF)
  colnames(beta_OR) <- c(var_individual, var_neighborhood)
  rownames(beta_OR) <- c("Odds Ratio", "CI_LB", "CI_UB", "bayes_p", "VIF")
  coeff_OR <- t(beta_OR) 
  coeff_OR <- coeff_OR[,c("Odds Ratio","bayes_p", "VIF")]
  # assign(ifelse(add_str == "", "orig_OR","reweight_OR"),coeff_OR)
  res <- coeff_OR
  res_noVIF <- coeff_OR[,-3]
  
  assign(ifelse(output_string == "PRETERM", "res_noVIF_PRE","res_noVIF_STILL"),res_noVIF)
}
xtable(cbind(res_noVIF_STILL[index,],res_noVIF_PRE[index,]))


#############################################################################
##################### Neighborhood cluster risk  (Table 3) ##################
#############################################################################

## See neigh_cluster_analysis.R for this table


#############################################################################
############# Detailed neighborhood cluster risk  (Table S2 and S3) #########
#############################################################################

## See neigh_cluster_analysis.R for this tables


#############################################################################
##################### Neighborhood comparison (Table S4) ####################
#############################################################################

logit <- function(x) log(x/(1-x))
var_neighborhood <- c("proportion Asian","proportion Hispanic","proportion Black", 
                      "proportion women", "poverty", 
                      "public assistance", "labor force", 
                      "recent birth", "high school grad", 
                      "college grad", "occupied housing", "housing violation", 
                      "violent crime", "nonviolent crime")
var_neighborhood_old <- c("Prop_Asian", "Prop_Hispanic_Latino", "Prop_Black", 
                          "Prop_women15to50", "Prop_Below100percPoverty_women15to50", 
                          "Prop_ReceivedPublicAssistanceIncome_women15to50", "Prop_InLaborForce_women16to50", 
                          "Prop_BirthsInPast12Months_women15to50", "Prop_HighschoolGrad_women15to50", 
                          "Prop_BachelorsDegree_women15to50", "log.occupied.housing", "log.violation", 
                          "log.violent", "log.nonviolent")
var_individual <- c("Hispanic","Black", "Asian", "multiple_birth", "age")
p1 <- length(var_individual)
p2 <- length(var_neighborhood)
col_to_be_scaled <- p1:(p1+p2)

philly_all_covariates$log.occupied.housing <- asinh(0.5*philly_all_covariates$Total_Occupied_Housing_Units)
philly_all_covariates$log.violation <- asinh(0.5*philly_all_covariates$violation)
philly_all_covariates$log.violent <- asinh(0.5*philly_all_covariates$violent)
philly_all_covariates$log.nonviolent <- asinh(0.5*philly_all_covariates$nonviolent)
# select only the covariates we use
philly_all_covariates <- (philly_all_covariates[,c(1,2,match(var_neighborhood_old, colnames(philly_all_covariates)))])

### we need this only for tractID_unique (tracts used in the analysis)
year <- 8; add_str = ""; add_str2 <- "_YEAR8"; output_string <- "STILLBIRTH"
tmp <-load(paste0(wdstr,"output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
tractID_unique <- tractID_unique_list[[year]] # used neighborhoods 
philly_all_covariates <- philly_all_covariates[which(philly_all_covariates$tractID %in% tractID_unique),] # select only the neighborhoods we use

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

################# get probabilities #################

year <- 8; add_str <- ""; add_str2 <- "_YEAR8"
# STILL
output_string <- "STILLBIRTH"
tmp <-load(paste0(wdstr,"phat_output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
phat_neig <- rowMeans(avg_phat_samples)  # already been thinned

data <- data.frame()

i <- 1
i1 <- IS1[i]; i2 <- IS2[i]
iB1 <- IBS1[i]; iB2 <- IBS2[i]

data <- rbind(data, philly_all_covariates[c(i1,i2),])

data$phat_ST[1:2] <- phat_neig[c(iB1,iB2)]
data$OR_ST[1:2] <- exp(logit(phat_neig[c(iB1,iB2)]) )
data$RelOR_ST[1:2] <- c(NA,exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(-1,1) ))) ## let's compute OR2/OR1

for(i in 2:4){
  i1 <- IS1[i]; i2 <- IS2[i]
  iB1 <- IBS1[i]; iB2 <- IBS2[i]
  data <- rbind(data, cbind(philly_all_covariates[c(i1,i2),], 
                            phat_ST = phat_neig[c(iB1,iB2)],
                            OR_ST = exp(logit(phat_neig[c(iB1,iB2)]) ),
                            RelOR_ST = c(NA, exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(-1,1) )))))
}

data$phat_PR <- NA
data$OR_PR <- NA
data$RelOR_PR <- NA

# PRE
output_string <- "PRETERM"
tmp <-load(paste0(wdstr,"phat_output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
phat_neig <- rowMeans(avg_phat_samples)  # already been thinned
for(i in 1:4){
  i1 <- IS1[i]; i2 <- IS2[i]
  iB1 <- IBS1[i]; iB2 <- IBS2[i]
  
  data$phat_PR[2*(i-1)+1:2] = phat_neig[c(iB1,iB2)]
  data$OR_PR[2*(i-1)+1:2] = exp(logit(phat_neig[c(iB1,iB2)]) )
  data$RelOR_PR[2*(i-1)+1:2] = c(NA, exp(sum( logit(phat_neig[c(iB1,iB2)]) * c(-1,1) )))
}

data_print <- data[,-c(1,2,18,21)] # do not print tractID, year and OR
colnames(data_print)[match(var_neighborhood_old,colnames(data_print))] <- var_neighborhood
data_print[, var_neighborhood] <- round(data_print[, var_neighborhood],2)
data_print[, c("phat_ST","RelOR_ST","phat_PR","RelOR_PR")] <- round(100*data_print[, c("phat_ST","RelOR_ST","phat_PR","RelOR_PR")],2)
for(col in c("phat_ST","RelOR_ST","phat_PR","RelOR_PR")){
  data_print[, col] <- paste0(data_print[, col],"%")
}
data_print[data_print == "NA%"] <- ""
data_print <- data_print[, c(15:18,1:14)]

View(t(data_print))
xtable(t(data_print))

#############################################################################
######## Predictive performance, CAR vs SMOTE vs reweight (table S5) ########
#############################################################################

method <- "CAR"; year <- 8; add_str2 <- "_YEAR8"
####### PRETERM
output_string <- "PRETERM"

add_str <- ""
tmp <-load(paste0(wdstr,"output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
noSMOTE_new <- c(auc[year], DIC[year], rmse[year], DIC_shifted[year], rmse_shifted[year], misclass[year], sens[year], spec[year])
add_str <- "a-reweight_"
tmp <-load(paste0(wdstr,"output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
reweigh_new <- c(auc[year], DIC[year], rmse[year], DIC_shifted[year], rmse_shifted[year], misclass[year], sens[year], spec[year])
add_str <- ""
tmp <- load(paste0(wdstr,"output_LOO_nogamma_rho_SMOTE_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
SMOTE_new <- c(auc[year], DIC[year], rmse[year], DIC_shifted[year], rmse_shifted[year], misclass[year], sens[year], spec[year])

res_pret_loo <- rbind(noSMOTE_new,SMOTE_new,reweigh_new)
colnames(res_pret_loo) <- c("AUC", "DIC","rmse","DIC*","rmse*","MR", "sens", "spec")
res_pret <- res_pret_loo[,-(2:5)]

####### STILLBIRTH
output_string <- "STILLBIRTH"

add_str <- ""
tmp <-load(paste0(wdstr,"output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
noSMOTE_new <- c(auc[year], DIC[year], rmse[year], DIC_shifted[year], rmse_shifted[year], misclass[year], sens[year], spec[year])
add_str <- "a-reweight_"
tmp <-load(paste0(wdstr,"output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
reweigh_new <- c(auc[year], DIC[year], rmse[year], DIC_shifted[year], rmse_shifted[year], misclass[year], sens[year], spec[year])
add_str <- ""
tmp <- load(paste0(wdstr,"output_LOO_nogamma_rho_SMOTE_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
SMOTE_new <- c(auc[year], DIC[year], rmse[year], DIC_shifted[year], rmse_shifted[year], misclass[year], sens[year], spec[year])

res_still_loo <- rbind(noSMOTE_new,SMOTE_new,reweigh_new)
colnames(res_still_loo) <- c("AUC", "DIC","rmse","DIC*","rmse*","MR", "sens", "spec")
res_still <- res_still_loo[,-(2:5)]

res <- t(rbind(res_still, res_pret))
xtable(res, digits = 3)


#############################################################################
######################## MCMC diagnostics (table S8) ########################
#############################################################################

var_neighborhood <- c("proportion Asian","proportion Hispanic","proportion Black",
                      "proportion women", "poverty", 
                      "public assistance", "labor force", 
                      "recent birth", "high school grad", 
                      "college grad", "occupied housing", "housing violation", 
                      "violent crime", "nonviolent crime")
var_individual <- c("Hispanic","Black", "Asian", "multiple birth", "age")
vars <- c(var_individual, var_neighborhood)

burnin <- 500
thin <- 5
mcmc_niter <- 1000*thin+burnin 
index_thinning <- seq(burnin, mcmc_niter, by = thin)

year <- 8; add_str = ""; add_str2 <- "_YEAR8"

for(output_string in c("PRETERM", "STILLBIRTH")){
  tmp <-load(paste0(wdstr,"output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
  
  cov_index <- 363 + 1:19
  
  betas <- cbind(tmp1$beta[,index_thinning, drop = FALSE],
                 tmp2$beta[,index_thinning, drop = FALSE])
  beta <- rowMeans(betas)
  
  MCMC_summary <- cbind(ES2[cov_index], beta, MCSE[cov_index])
  
  tmp1 <- output_list1[[year]]; tmp2 <- output_list2[[year]]
  mat1 <- as.matrix(t(rbind(tmp1$alpha, tmp1$beta, tmp1$rho, 
                            tmp1$a0, tmp1$b0,
                            tmp1$tau2_alpha, tmp1$tau2_beta)))
  mat2 <- as.matrix(t(rbind(tmp2$alpha, tmp2$beta, tmp2$rho,
                            tmp2$a0, tmp2$b0,
                            tmp2$tau2_alpha, tmp2$tau2_beta)))
  mat_cov <- rbind(mat1,mat2)[,cov_index]
  
  assign(paste0("MCMC_summary_",output_string),MCMC_summary)
  assign(paste0("mat_cov_",output_string),mat_cov)
  
  # assign(paste0("tmp1_",output_string),tmp1)
  # assign(paste0("tmp2_",output_string),tmp2)
}

MCMC_summary <- cbind(MCMC_summary_STILLBIRTH,MCMC_summary_PRETERM)
colnames(MCMC_summary) <- c("ESS", "log-odds", "MCSE", "ESS", "log-odds", "MCSE")
rownames(MCMC_summary) <- vars

xtable(MCMC_summary, digits  = c(NA,1,3,4,1,3,4))

