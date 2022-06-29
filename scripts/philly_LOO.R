## please call this script from meta_script.R (or assign the values initialized in meta_script.R before running this).

index_thinning <- seq(burnin, mcmc_niter, by = thin)

if(is.factor(data.preg.tractcov.atleast10$YEAR)){
  years = levels(data.preg.tractcov.atleast10$YEAR)
} else {
  years = sort(unique(data.preg.tractcov.atleast10$YEAR))
}
n_years = length(years)
n = nrow(data.preg.tractcov.atleast10)

num.reps = n_years


# Vectors to store the results
rmse = rep(NA,num.reps)       # Holds the overall RMSE.
misclass = rep(NA,num.reps)   # Holds the overall misclassification rate.
sens = rep(NA,num.reps)       # Holds the sensitivity (true positive rate)
spec = rep(NA,num.reps)       # Holds the specificity (true negative rate)
auc = rep(NA,num.reps)        # Holds the area under the ROC curve
DIC = rep(NA,num.reps)        # Holds the overall DIC
WAIC1 = rep(NA,num.reps)        # Holds the overall WAIC
WAIC2 = rep(NA,num.reps)        # Holds the overall WAIC
thresholds = rep(NA,num.reps) # Holds the thresholds used for sens and spec
medianrho = rep(NA,num.reps)  # Holds the median rho found from the new MCMC
rmse_shifted = rep(NA,num.reps)
DIC_shifted = rep(NA,num.reps)
WAIC1_shifted = rep(NA,num.reps)
WAIC2_shifted = rep(NA,num.reps)
prauc <- numeric(num.reps)
prop_cases_train <- numeric(num.reps)
prop_cases_test <- numeric(num.reps)

output_list1 <- list()
output_list2 <- list()
roc_output_list <- list()
XtX_list <- list()
my_scale_list <- list()
scale_means_list <- list()
scale_sds_list <- list()
tractID_unique_list <- list()
alpha_reweight_list <- list()
beta_reweight_list <- list()

# new phat's
avg_phat_samples_list <- list()
avg_yhat_samples_list <- list()
avg_phat2_samples_list <- list()
avg_phat2shift_samples_list <- list()
avg_phat_shift_samples_list <- list()
# old phat's
mean_phat_is_list <- list()
mean_phat2_is_list <- list()
mean_phat2shift_is_list <- list()
mean_phat_all_list <- list()
mean_phat2_all_list <- list()
mean_phat2shift_all_list <- list()
phat_neig_list <- list()
phat_neigshift_list <- list()

gamma_opts <- array(NA, dim = c(num.reps, 3))
f_opts <- array(NA, dim = c(num.reps, 3))
gamma_opts_new <- array(NA, dim = c(num.reps, 2))
f_opts_new <- array(NA, dim = c(num.reps, 2))

### We do not include WHITE or Prop_White to avoid collinearity issues
var_neighborhood <- c("Prop_Asian", "Prop_Hispanic_Latino", "Prop_Black", 
                      "Prop_women15to50", "Prop_Below100percPoverty_women15to50", 
                      "Prop_ReceivedPublicAssistanceIncome_women15to50", "Prop_InLaborForce_women16to50", 
                      "Prop_BirthsInPast12Months_women15to50", "Prop_HighschoolGrad_women15to50", 
                      "Prop_BachelorsDegree_women15to50", "log.occupied.housing", "log.violation", 
                      "log.violent", "log.nonviolent")
var_individual <- c("HISPANIC", "NH_AFAM", "NH_ASIAN", "MULTIPLE_BIRTH", "ENC_AGE_DEC")
var_individual_SMOTE <- c("RACE", "MULTIPLE_BIRTH", "ENC_AGE_DEC")
if(interaction_bool == TRUE){
  var_individual <- c(var_individual, "black_interaction", "white_interaction", "asian_interaction", "hispanic_interaction")
}

p1 <- length(var_individual)
p2 <- length(var_neighborhood)
ps_01 <- var_individual[-p1]

if(neigh_bool == TRUE){
  p_used <- p1+p2
  var_used <- c(var_individual, var_neighborhood)
  var_used_SMOTE <- c(var_individual_SMOTE, var_neighborhood)
  col_to_be_scaled_names <- c("ENC_AGE_DEC",var_neighborhood) # of the individual, we only scale AGE
} else {
  p_used <- p1
  var_used <- var_individual
  var_used_SMOTE <- var_individual_SMOTE
  col_to_be_scaled_names <- "ENC_AGE_DEC"                     # of the individual, we only scale AGE
}
sig.cov = array(NA,dim = c(p_used, num.reps))
rownames(sig.cov) <- var_used

get_covmat_scaled <- function(data.preg.tractcov.atleast10TMP, ps_01, col_to_be_scaled_names, var_used){
  covmat <- data.matrix(data.preg.tractcov.atleast10TMP[, var_used])
  if(all(unique(covmat[,ps_01]) > 0)){
    covmat <- covmat - 1 # because factors were turned into 1/2 and we want 0/1
  }
  covmat_scaled <- covmat
  for(k in col_to_be_scaled_names){
    covmat_scaled[,k] <- my_scale(covmat[,k],k)
  }
  return(covmat_scaled)
}

for(j in 1:ncol(data.preg.tractcov.atleast10)){
  if(colnames(data.preg.tractcov.atleast10)[j] != "tractID"){
    data.preg.tractcov.atleast10[,j] <- as.numeric(as.character(data.preg.tractcov.atleast10[,j]))
  }
}

set.seed(123456)
if(newseed_bool){
  set.seed(12345678)
}
seed_smote <- sample(x = 1e+6, size = 1)
seeds_aftersmote <- sample(x = 1e+6, size = num.reps)
for(i_year in c(num.reps, 1:(num.reps-1))){
  ## find test and training indices
  test.indices = which(data.preg.tractcov.atleast10$YEAR == years[i_year])
  training.indices = setdiff(c(1:n),test.indices) # the remaining ones
  
  ## create training set
  data.preg.tractcov.atleast10TMP = data.preg.tractcov.atleast10[training.indices,]
  
  if(sum(is.na(data.preg.tractcov.atleast10TMP))>0){
    warning("Some NA present in training set, iter", i_year)
    indexNA <- which(apply(data.preg.tractcov.atleast10TMP[,var_used], 
                           MARGIN = 1, function(x) any(is.na(x))))
    data.preg.tractcov.atleast10TMP <- data.preg.tractcov.atleast10TMP[-indexNA,]
  }
  
  overall_prob_train <- mean(data.preg.tractcov.atleast10TMP[,outcome])
  data.preg.tractcov.atleast10TMP <- binary_to_factor(data.preg.tractcov.atleast10TMP)

  if(smote_bool == TRUE){
    data.preg.tractcov.atleast10TMP[,outcome] <- as.factor(data.preg.tractcov.atleast10TMP[,outcome])
    data.preg.tractcov.atleast10TMP[,"YEAR"] <- as.factor(data.preg.tractcov.atleast10TMP[,"YEAR"])
    
    over <- 200
    under <- 200
    
    set.seed(seed_smote)
    if(outcome == "PRETERM"){
      data_smoted <- performanceEstimation::smote(PRETERM ~ .,
                                                  data = data.preg.tractcov.atleast10TMP[,c("LONGITUDE", "LATITUDE","YEAR", var_individual_SMOTE,outcome)],
                                                  perc.over = over/100, perc.under = under/100)
    } else if(outcome == "STILLBIRTH"){
      data_smoted <- performanceEstimation::smote(STILLBIRTH ~ .,
                                                  data = data.preg.tractcov.atleast10TMP[,c("LONGITUDE", "LATITUDE","YEAR", var_individual_SMOTE,outcome)],
                                                  perc.over = over/100, perc.under = under/100)
    } else if(outcome == "CSECTION"){
      data_smoted <- performanceEstimation::smote(CSECTION ~ .,
                                                  data = data.preg.tractcov.atleast10TMP[,c("LONGITUDE", "LATITUDE","YEAR", var_individual_SMOTE,outcome)],
                                                  perc.over = over/100, perc.under = under/100)
    } else {
      warning("Outcome ", outcome," not supported.")
    }
    data.preg.tractcov.atleast10TMP[,"YEAR"] <- as.numeric(as.character(data.preg.tractcov.atleast10TMP[,"YEAR"]))
    data_smoted[,"YEAR"] <- as.numeric(as.character(data_smoted[,"YEAR"]))
    set.seed(seeds_aftersmote[i_year])
    
    ### need to deal with the geographic location and tractID
    pts <- SpatialPoints(data_smoted[ ,c('LONGITUDE', 'LATITUDE')])
    databytracts <- sp::over(tracts, pts, returnList = T)
    datavector <- unlist(databytracts)
    datatractlength <- lapply(databytracts, length)
    dataind.temp <- c(0, cumsum(datatractlength))
    levels_region <- levels(tracts@data$GEOID10)[tracts@data$GEOID10]
    mapping <- data.frame(tract = NA, tractID = NA, row_index = datavector)
    for(i in 1:length(databytracts)){
      if(datatractlength[i]>0){
        temp.index <- (dataind.temp[i] + 1): (dataind.temp[i + 1])
        mapping$tract[temp.index] = i
        mapping$tractID[temp.index] = levels_region[i]
      }
    }
    
    
    ## we should exclude the rows that are NOT in the tracts that we are considering
    to_keep <- which(mapping$tractID %in% tractID_tokeep)
    to_throw <- which(!(mapping$tractID %in% tractID_tokeep))
    mapping_keep <- mapping[to_keep,]
    data_smoted_keep <- data_smoted[mapping_keep$row_index,] 
    ## now data_smoted_keep is sorted according to mapping (tract and tractID)
    ## now mapping$row_index does not correspond to the row index of data_smoted_keep, but of the old data_smoted
    data_smoted_keep <- cbind(data_smoted_keep, mapping_keep$tractID)
    colnames(data_smoted_keep)[ncol(data_smoted_keep)] <- "tractID"
    

    # let's create data.preg.tractcov.atleast10TMP
    if(neigh_bool == TRUE){
      ## let's create data.preg.tractcov.atleast10TMP by matching the neighborhood covariates 
      index <- match(interaction(data_smoted_keep[,c("tractID","YEAR")]), interaction(data.preg.tractcov.atleast10TMP[,c("tractID","YEAR")]))
      data.preg.tractcov.atleast10TMP2 <- cbind(data_smoted_keep, data.preg.tractcov.atleast10TMP[index, var_neighborhood])
    } else {
      data.preg.tractcov.atleast10TMP2 <- data_smoted_keep
    }
    if(sum(is.na(data.preg.tractcov.atleast10TMP2))>0){
      ## some NA can be created if there was no tractID-YEAR combination in the original set.
      # warning("Still NA present - removing them in iter",i_year)
      indexNA <- which(apply(data.preg.tractcov.atleast10TMP2[,var_used_SMOTE], # maybe you could just use "var_used" here
                             MARGIN = 1, function(x) any(is.na(x))))
      data.preg.tractcov.atleast10TMP2 <- data.preg.tractcov.atleast10TMP2[-indexNA,]
      mapping_keep <- mapping_keep[-indexNA,]
    }


    # if some neighborhood are missing, we add them back from the original dataset
    if(length(unique(mapping_keep$tractID)) < length(tractID_tokeep)){
      tractID_kept <- unique(mapping_keep$tractID)
      tractID_missing <- tractID_tokeep[which(!(tractID_tokeep %in% tractID_kept))]
      index_reuse <- which(data.preg.tractcov.atleast10TMP$tractID %in% tractID_missing)
      if(length(tractID_missing)>0 & length(index_reuse)==0){
        warning("data.preg.tractcov.atleast10TMP does not contain data from tractID_missing.")
      }
      ## we will just add the data that was excluded (probably from the majority class)
      add <- data.preg.tractcov.atleast10TMP[index_reuse,c("LONGITUDE", "LATITUDE","YEAR", var_individual_SMOTE,outcome,"tractID",var_neighborhood)]
      data.preg.tractcov.atleast10TMP2 <- rbind(data.preg.tractcov.atleast10TMP2, add)
      mapping_add <- data.frame(tract = NA, 
                                tractID = data.preg.tractcov.atleast10TMP[index_reuse, "tractID"], 
                                row_index = (nrow(data.preg.tractcov.atleast10TMP2)-nrow(add)+1):nrow(data.preg.tractcov.atleast10TMP2))
      for(i in 1:length(index_reuse)){
        mapping_add[i,"tract"] <- which(mapping_add[i,"tractID"] == levels_region)
      }
      mapping_keep <- rbind(mapping_keep, mapping_add)
    }
    data.preg.tractcov.atleast10TMP <- data.preg.tractcov.atleast10TMP2
    data.preg.tractcov.atleast10TMP[,outcome] = as.numeric(as.character(data.preg.tractcov.atleast10TMP[,outcome]))
    
    ## now let's create new_mapping. mapping is sorted like data.preg.tractcov.atleast10TMP
    mapping <- mapping_keep 
    # columns are tract, tractID, row_index
    tractID_unique <- unique(mapping$tractID) # these are not sorted, they're in the order of appearance
    tractID_unique_list[[i_year]] <- tractID_unique
    mapping$new_tract <- match(mapping$tractID, tractID_unique) 
    # NOW columns are tract, tractID, row_index, new_tract
    new_mapping <- mapping[,c(4,2)]
    names(new_mapping)[1] <- "tract"
    if(noRE_bool){
      ## this is only for single intercept:
      new_mapping$tract <- 1 
    }
    mapping <- new_mapping
    # NOW columns are tract(new_tract), tractID
    
  } else {
    unik_tractID <- unique(data.preg.tractcov.atleast10TMP$tractID)
    unik_yearID <- unique(data.preg.tractcov.atleast10TMP$YEAR)
    mapping <- data.frame(tract = NA, year = NA, 
                          tractID = data.preg.tractcov.atleast10TMP$tractID, 
                          yearID = data.preg.tractcov.atleast10TMP$YEAR)
    mapping$tract = match(mapping$tractID, unik_tractID)
    mapping$year = match(mapping$yearID, unik_yearID)
    tractID_unique <- unique(mapping$tractID)
    tractID_unique_list[[i_year]] <- tractID_unique
    if(noRE_bool){
      ## this is only for single intercept:
      mapping$tract <- 1 
    }
  }
  
  data.preg.tractcov.atleast10TMP <- factor_to_binary(data.preg.tractcov.atleast10TMP)
  data.preg.tractcov.atleast10TMP <- data.preg.tractcov.atleast10TMP[,c("tractID",outcome,var_used)] # in this way it's the right order (and excludes WHITE)
  nalpha = nrow(w.sym)
  
  if(length(unique(mapping$tractID)) < length(tractID_tokeep)){
    warning("No datapoints for some neighborhoods in iter ",i_year,"- skipping this one.")
    next
  }
  if( (nalpha > 1) & (nalpha != length(unique(mapping$tractID))) ){
    warning("Number of neighborhoods in w.sym is different from the one in the data for iter ",i_year,"- skipping this one.")
    next
  }
  
  ## create covariate matrix
  covmat <- data.matrix(data.preg.tractcov.atleast10TMP[, var_used])
  for(v in ps_01){
    if(all(unique(covmat[,v]) > 0)){ 
      covmat[,v] <- covmat[,v] - 1 # because factors were turned into 1/2 and we want 0/1
    }
  }
  covmat_scaled <- covmat
  covmat_scaled[,col_to_be_scaled_names] <- scale(covmat[,col_to_be_scaled_names])
  scale_means <- apply(X = covmat, MARGIN = 2, mean, na.rm = T); names(scale_means) = var_used
  scale_sds <- apply(X = covmat, MARGIN = 2, sd, na.rm = T); names(scale_sds) = var_used
  my_scale <- function(x, k){
    (x - scale_means[k])/scale_sds[k]
  }
  my_scale_list[[i_year]] <- my_scale
  scale_means_list[[i_year]] <- scale_means
  scale_sds_list[[i_year]] <- scale_sds

  prop_cases_train[i_year] <- mean(data.preg.tractcov.atleast10TMP[,outcome])
  
  w.sym_temp <-  w.sym
  if(prior == "CAR"){
    rho <- 0.9        # this is a starting point
    sample_rho <- TRUE
  } else if(prior == "independent"){
    rho <- 0
    sample_rho <- FALSE
  }
  if(noRE_bool){
    w.sym_temp <- matrix(1, nrow = 1, ncol = 1)
    rho <- 0
    sample_rho <- FALSE
  }

  ## Construct the weights for the training set
  n.train = nrow(data.preg.tractcov.atleast10TMP)
  num.cases = sum(data.preg.tractcov.atleast10TMP[,outcome])
  num.controls = n.train-num.cases
  ratio = num.cases/num.controls
  if(reweight_bool){
    if(alpha_vs_beta_bool){
      alpha_reweight = ratio
      beta_reweight = 1
    } else {
      alpha_reweight = 1
      beta_reweight = 1/ratio
    }
  } else {
    alpha_reweight = 1
    beta_reweight = 1
  }
  alpha_reweight_list[[i_year]] <- alpha_reweight
  beta_reweight_list[[i_year]] <- beta_reweight
  
  ptm <- proc.time()
  tmp1 <- mcmc(y = data.preg.tractcov.atleast10TMP[,outcome], 
            mapping = mapping,
            X = covmat_scaled,
            niter = mcmc_niter,
            w.sym = w.sym_temp, rho = rho, a_rho = 1, b_rho = 1,
            tau2_alpha = 0.1, tau2_beta = 0.1,
            alpha_reweight = alpha_reweight,
            beta_reweight = beta_reweight,
            sample_rho = sample_rho)
  tmp2 <- mcmc(y = data.preg.tractcov.atleast10TMP[,outcome], 
               mapping = mapping,
               X = covmat_scaled,
               niter = mcmc_niter,
               w.sym = w.sym_temp, rho = rho, a_rho = 1, b_rho = 1,
               tau2_alpha = 0.1, tau2_beta = 0.1,
               alpha_reweight = alpha_reweight,
               beta_reweight = beta_reweight,
               sample_rho = sample_rho)
  ptm <- proc.time() - ptm
  output_list1[[i_year]] <- tmp1
  output_list2[[i_year]] <- tmp2
  
  mat1 <- as.matrix(t(rbind(tmp1$alpha, tmp1$beta, tmp1$rho, 
                            tmp1$a0, tmp1$b0,
                            tmp1$tau2_alpha, tmp1$tau2_beta)))
  mat2 <- as.matrix(t(rbind(tmp2$alpha, tmp2$beta, tmp2$rho,
                            tmp2$a0, tmp2$b0,
                            tmp2$tau2_alpha, tmp2$tau2_beta)))
  ES2 <- apply( rbind(mat1,mat2), MARGIN = 2, FUN = LaplacesDemon::ESS)
  MCSE <- apply( rbind(mat1,mat2), MARGIN = 2, FUN = LaplacesDemon::MCSE)
    
  mat1_thin <- mat1[index_thinning,, drop = FALSE]
  mat2_thin <- mat2[index_thinning,, drop = FALSE]
  ES2_thin <- apply( rbind(mat1_thin,mat2_thin), MARGIN = 2, FUN = LaplacesDemon::ESS)
  MCSE_thin <- apply( rbind(mat1_thin,mat2_thin), MARGIN = 2, FUN = LaplacesDemon::MCSE)  
  
  tmp <- list(alpha = cbind(tmp1$alpha[,index_thinning, drop = FALSE],
                          tmp2$alpha[,index_thinning, drop = FALSE]),
            beta = cbind(tmp1$beta[,index_thinning, drop = FALSE],
                         tmp2$beta[,index_thinning, drop = FALSE]),
            rho = c(tmp1$rho[index_thinning],tmp2$rho[index_thinning]))

  alpha_pm <- rowMeans(tmp$alpha) # posterior mean
  beta_pm <- rowMeans(tmp$beta) # posterior mean
  beta_est <- beta_pm
  medianrho[i_year] <- median(tmp$rho)
  
  beta_CI <- apply(tmp$beta, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975)) 
  CIprod <- apply(beta_CI, MARGIN = 2, FUN = prod)
  contains_zero <- ifelse(CIprod < 0, 1, 0)
  sig.cov[,i_year] <- 1-contains_zero
  
  ## check collinearity
  X = covmat_scaled
  XtX <- t(X) %*% X
  XtX_list[[i_year]] <- XtX
  
  ############################################################
  ### all sample - compute probabilities for neighborhood analysis
  data.preg.tractcov.atleast10TMP <- data.preg.tractcov.atleast10
  covmat_scaled <- get_covmat_scaled(data.preg.tractcov.atleast10TMP, ps_01, col_to_be_scaled_names, var_used)
  temp = match(data.preg.tractcov.atleast10TMP$tractID, mapping$tractID)
  index_mapping <- mapping$tract[temp]

  ### computing phat using posterior samples (used currently)
  re_samples <- tmp$alpha[index_mapping,]
  yhat_samples <- covmat_scaled %*% as.matrix(tmp$beta) + re_samples
  phat_samples <- exp(yhat_samples) / (1 + exp(yhat_samples))

  avg_phat_samples <- array(NA, dim = c(length(tractID_unique),n_samples)) # this takes the expectation of the phats
  weights <- numeric(length(tractID_unique))
  avg_yhat_samples <- array(NA, dim = c(length(tractID_unique),n_samples))
  for(i in 1:length(tractID_unique)){
    index <- which(data.preg.tractcov.atleast10TMP$tractID == tractID_unique[i])
    
    avg_phat_samples[i,] <- colMeans(phat_samples[index,,drop=FALSE])
    weights[i] <- length(index)/nrow(data.preg.tractcov.atleast10TMP)
    avg_yhat_samples[i,] <- colMeans(yhat_samples[index,,drop=FALSE])
  }
  avg_phat2_samples <- exp(avg_yhat_samples) / (1 + exp(avg_yhat_samples)) # this takes the expectation of the yhat and then transforms

  # find gammas for recalibration
  myfun1 <- function(gamma){
    tmp <- mean( exp(yhat_samples - gamma) / (1 + exp(yhat_samples - gamma))  )
    (tmp-overall_prob_train)^2
  }
  myfun2 <- function(gamma){
    tmp <- sum( weights * rowMeans(exp(avg_yhat_samples - gamma) / (1 + exp(avg_yhat_samples - gamma))) )
    (tmp-overall_prob_train)^2
  }
  # for avg_phat_samples (avg_phat_shift_samples)
  tmp_opt1 <- optimize(myfun1, c(-10,10))
  gamma_opt1 <- tmp_opt1$minimum
  f_opt1 <- tmp_opt1$objective
  # for avg_phat2_samples (avg_phat2shift_samples)
  tmp_opt2 <- optimize(myfun2, c(-10,10))
  gamma_opt2 <- tmp_opt2$minimum
  f_opt2 <- tmp_opt2$objective

  gamma_opts_new[i_year,1] <- gamma_opt1
  f_opts_new[i_year,1] <- f_opt1
  gamma_opts_new[i_year,2] <- gamma_opt2
  f_opts_new[i_year,2] <- f_opt2

  avg_phat_shift_samples <- array(NA, dim = c(length(tractID_unique),n_samples))
  for(i in 1:length(tractID_unique)){
    index <- which(data.preg.tractcov.atleast10TMP$tractID == tractID_unique[i])
    temp = exp(yhat_samples-gamma_opt1) / (1 + exp(yhat_samples-gamma_opt1))
    avg_phat_shift_samples[i,] <- colMeans( temp[index,,drop=FALSE])
  }
  avg_phat2shift_samples <- exp(avg_yhat_samples - gamma_opt2) / (1 + exp(avg_yhat_samples - gamma_opt2))

  avg_phat_samples_list[[i_year]] <- avg_phat_samples
  avg_yhat_samples_list[[i_year]] <- avg_yhat_samples
  avg_phat2_samples_list[[i_year]] <- avg_phat2_samples
  avg_phat_shift_samples_list[[i_year]] <- avg_phat_shift_samples
  avg_phat2shift_samples_list[[i_year]] <- avg_phat2shift_samples

  ### computing phat after computing posterior mean (discarded)
  re_est = alpha_pm[mapping$tract[temp]] 
  yhat_all <- covmat_scaled %*% matrix(beta_est, ncol = 1) + re_est
  phat_all <- exp(yhat_all) / (1 + exp(yhat_all))
    
  mean_phat_all <- numeric(length(tractID_unique)) # this takes the expectation of the phats
  mean_phat2_all <- numeric(length(tractID_unique)) # this takes the expectation of the yhat and then transforms
  mean_phat2shift_all <- numeric(length(tractID_unique))
    
  weights <- numeric(length(tractID_unique))
  avg_yhat_all <- numeric(length(tractID_unique))
  for(i in 1:length(tractID_unique)){
    index <- which(data.preg.tractcov.atleast10TMP$tractID == tractID_unique[i])
    
    mean_phat_all[i] <- mean(phat_all[index])
    weights[i] <- length(index)/nrow(data.preg.tractcov.atleast10TMP)
    avg_yhat_all[i] <- mean(yhat_all[index])
  }
  myfun <- function(gamma){
    tmp <- sum(weights * exp(avg_yhat_all - gamma) / (1 + exp(avg_yhat_all - gamma)))
    (tmp-overall_prob_train)^2
  }
  tmp_opt <- optimize(myfun, c(-10,10))
  gamma_opt <- tmp_opt$minimum
  f_opt <- tmp_opt$objective
  gamma_opts[i_year,2] <- gamma_opt
  f_opts[i_year,2] <- f_opt
    
  mean_phat2_all <- exp(avg_yhat_all) / (1 + exp(avg_yhat_all))
  mean_phat2shift_all <- exp(avg_yhat_all - gamma_opt) / (1 + exp(avg_yhat_all - gamma_opt))

  mean_phat_all_list[[i_year]] <- mean_phat_all
  mean_phat2_all_list[[i_year]] <- mean_phat2_all
  mean_phat2shift_all_list[[i_year]] <- mean_phat2shift_all
    
  ### Create test set for predictive performance
  data.preg.tractcov.atleast10TMP = data.preg.tractcov.atleast10[test.indices,]
  prop_cases_test[i_year] <- mean(data.preg.tractcov.atleast10TMP[,outcome])
  covmat_scaled <- get_covmat_scaled(data.preg.tractcov.atleast10TMP, ps_01, col_to_be_scaled_names, var_used)
  # some tracts will have 0 pregnancies out of sample, so they won't appear in the new 
  temp = match(data.preg.tractcov.atleast10TMP$tractID, mapping$tractID)

  re_est = alpha_pm[mapping$tract[temp]] # this is equivalent to reordering mapping so that the order matches the out of sample
  yhat_os <- covmat_scaled %*% matrix(beta_est, ncol = 1) + re_est
  phat_os <- exp(yhat_os) / (1 + exp(yhat_os))
  phat_os_shifted <- exp(yhat_os - gamma_opt) / (1 + exp(yhat_os - gamma_opt))

  y_os <- data.preg.tractcov.atleast10TMP[, outcome]
  prauc[i_year] = PRAUC(y_pred = as.numeric(phat_os), y_true = as.numeric(y_os))
  rmse[i_year] <- mean((y_os - phat_os)^2)
  rmse_shifted[i_year] <- mean((y_os - phat_os_shifted)^2)
  f1 = roc(as.numeric(y_os) ~ as.numeric(phat_os), levels = c(0,1), direction = "<")      # Obtain the AUC (area under the curve)
  roc_output_list[[i_year]] <- f1
  auc[i_year] = as.numeric(auc(f1))
  threshold_Youden <- f1$thresholds[which.max(f1$sensitivities + f1$specificities - 1)]
  thresholds[i_year] = threshold_Youden
  misclass[i_year] <- mean(ifelse(y_os - phat_os >= threshold_Youden, 1, 0))
  correct.pos <- which(phat_os >= threshold_Youden & y_os==1)         # Sensitivity (true positive rate)
  sens[i_year] <- length(correct.pos)/length(which(y_os==1))
  correct.neg <- which(phat_os < threshold_Youden & y_os==0)          # Specificity (true negative rate)
  spec[i_year] <- length(correct.neg)/length(which(y_os==0))
    
  # compute DIC (this is computed on out of sample data)
  alphai_samples <- tmp$alpha[mapping$tract[temp],, drop = FALSE] + covmat_scaled %*% tmp$beta
  pi_samples <- exp(alphai_samples)/(1+exp(alphai_samples))
  pi_pm <- rowMeans(pi_samples)
  yi <- data.preg.tractcov.atleast10TMP[,outcome]
  logp_bayes <- sum(yi * log(pi_pm) + (1-yi) * log(1 - pi_pm))
  post_logp <- sum(yi * rowMeans(log(pi_samples)) + (1-yi) * rowMeans(log(1-pi_samples)))
  p_DIC <- 2 * (logp_bayes - post_logp)
  DIC[i_year] <- -2*logp_bayes + 2*p_DIC

  # compute WAIC (two forms, we discarded WAIC1)
  LPPD <- sum(log(rowMeans( yi * pi_samples + (1-yi) * (1 - pi_samples) )))
  p_WAIC1 <- 2 * sum(  
    log(rowMeans( yi * pi_samples + (1-yi) * (1 - pi_samples) )) - 
      rowMeans( yi * log(pi_samples) + (1-yi) * log(1 - pi_samples) )
  )
  p_WAIC2 <- sum(apply(yi * log(pi_samples) + (1-yi) * log(1-pi_samples), 1, var))
  WAIC1[i_year] <- -2*(LPPD - p_WAIC1)
  WAIC2[i_year] <- -2*(LPPD - p_WAIC2)
    
  # same but using shifted (calibrated) probabilities
  pi_samples_shifted <- exp(alphai_samples - gamma_opt)/(1+exp(alphai_samples - gamma_opt))
  pi_pm_shifted <- rowMeans(pi_samples_shifted)
  logp_bayes_shifted <- sum(yi * log(pi_pm_shifted) + (1-yi) * log(1 - pi_pm_shifted))
  post_logp_shifted <- sum(yi * rowMeans(log(pi_samples_shifted)) + (1-yi) * rowMeans(log(1-pi_samples_shifted)))
  p_DIC_shifted <- 2 * (logp_bayes_shifted - post_logp_shifted)
  DIC_shifted[i_year] <- -2*logp_bayes_shifted + 2*p_DIC_shifted

  LPPD_shifted <- sum(log(rowMeans( yi * pi_samples_shifted + (1-yi) * (1 - pi_samples_shifted) )))
  p_WAIC1_shifted <- 2 * sum(  
    log(rowMeans( yi * pi_samples_shifted + (1-yi) * (1 - pi_samples_shifted) )) - 
      rowMeans( yi * log(pi_samples_shifted) + (1-yi) * log(1 - pi_samples_shifted) )
  )
  p_WAIC2_shifted <- sum(apply(yi * log(pi_samples_shifted) + (1-yi) * log(1-pi_samples_shifted), 1, var))
  WAIC1_shifted[i_year] <- -2*(LPPD_shifted - p_WAIC1_shifted)
  WAIC2_shifted[i_year] <- -2*(LPPD_shifted - p_WAIC2_shifted)
    
  ### in sample - compute probabilities for neighborhood analysis
  data.preg.tractcov.atleast10TMP <- data.preg.tractcov.atleast10[training.indices,]
  covmat_scaled <- get_covmat_scaled(data.preg.tractcov.atleast10TMP, ps_01, col_to_be_scaled_names, var_used)
  temp = match(data.preg.tractcov.atleast10TMP$tractID, mapping$tractID)
  re_est = alpha_pm[mapping$tract[temp]]
  yhat_is <- covmat_scaled %*% matrix(beta_est, ncol = 1) + re_est
  phat_is <- exp(yhat_is) / (1 + exp(yhat_is))

  mean_phat_is <- numeric(length(tractID_unique))
  mean_phat2_is <- numeric(length(tractID_unique))
  mean_phat2shift_is <- numeric(length(tractID_unique))

  weights <- numeric(length(tractID_unique))
  avg_yhat_is <- numeric(length(tractID_unique))
  for(i in 1:length(tractID_unique)){
    index <- which(data.preg.tractcov.atleast10TMP$tractID == tractID_unique[i])
    mean_phat_is[i] <- mean(phat_is[index])
    weights[i] <- length(index)/nrow(data.preg.tractcov.atleast10TMP)
    avg_yhat_is[i] <- mean(yhat_is[index])
  }
  myfun <- function(gamma){
    tmp <- sum(weights * exp(avg_yhat_is - gamma) / (1 + exp(avg_yhat_is - gamma)))
    (tmp-overall_prob_train)^2
  }
  tmp_opt <- optimize(myfun, c(-10,10))
  gamma_opt <- tmp_opt$minimum
  f_opt <- tmp_opt$objective
  gamma_opts[i_year,1] <- gamma_opt
  f_opts[i_year,1] <- f_opt

  mean_phat2_is <- exp(avg_yhat_is) / (1 + exp(avg_yhat_is))
  mean_phat2shift_is <- exp(avg_yhat_is - gamma_opt) / (1 + exp(avg_yhat_is - gamma_opt))

  mean_phat_is_list[[i_year]] <- mean_phat_is
  mean_phat2_is_list[[i_year]] <- mean_phat2_is
  mean_phat2shift_is_list[[i_year]] <- mean_phat2shift_is
  
  ### all sample - ignoring individual covariates
  if(neigh_bool == TRUE){
    data.preg.tractcov.atleast10TMP <- data.preg.tractcov.atleast10
    covmat_scaled <- get_covmat_scaled(data.preg.tractcov.atleast10TMP, ps_01, col_to_be_scaled_names, var_used)
    covmat_scaled <- covmat_scaled[,p1 + 1:p2] # let's keep only var_neighborhood
    temp = match(tractID_unique, data.preg.tractcov.atleast10TMP$tractID)
    temp2 = match(tractID_unique, mapping$tractID)
    re_est = alpha_pm[mapping$tract[temp2]] 
    yhat_neig <- covmat_scaled[temp,] %*% matrix(beta_est[p1 + 1:p2], ncol = 1) + re_est
    
    weights <- numeric(length(tractID_unique))
    for(i in 1:length(tractID_unique)){
      weights[i] <- mean(data.preg.tractcov.atleast10TMP$tractID == tractID_unique[i])
    }
    
    myfun <- function(gamma){
      tmp <- sum(weights * exp(yhat_neig - gamma) / (1 + exp(yhat_neig - gamma)))
      (tmp-overall_prob_train)^2
    }
    tmp_opt <- optimize(myfun, c(-10,10))
    gamma_opt <- tmp_opt$minimum
    f_opt <- tmp_opt$objective
    gamma_opts[i_year,3] <- gamma_opt
    f_opts[i_year,3] <- f_opt
    
    phat_neig <- exp(yhat_neig) / (1 + exp(yhat_neig))
    phat_neigshift <- exp(yhat_neig - gamma_opt) / (1 + exp(yhat_neig - gamma_opt))
    
    phat_neig_list[[i_year]] <- phat_neig
    phat_neigshift_list[[i_year]] <- phat_neigshift
  }

}

varlist <- c("outcome","prior","interaction_bool","smote_bool","neigh_bool", "var_used", "p_used",
             "misclass", "sens", "spec", "auc", "sig.cov", "DIC", "rmse", "thresholds",
             "DIC_shifted", "rmse_shifted","index_thinning",
             "output_list1","output_list2","XtX_list","my_scale_list", "scale_means_list", "scale_sds_list",
             "tractID_unique_list",
             # new phat's
             "gamma_opts_new","f_opts_new",
             "avg_phat_samples_list","avg_yhat_samples_list","avg_phat2_samples_list",
             "avg_phat2shift_samples_list","avg_phat_shift_samples_list"
             # old phat's
             "mean_phat_is_list","mean_phat_all_list", "phat_neig_list",
             "mean_phat2_is_list","mean_phat2_all_list", 
             "mean_phat2shift_is_list","mean_phat2shift_all_list", "phat_neigshift_list",
             #
             "alpha_reweight_list","beta_reweight_list",
             "roc_output_list",
             "WAIC1_shifted", "WAIC2_shifted", "WAIC1", "WAIC2",
             "medianrho", "gamma_opts", "f_opts",
             "prauc","prop_cases_test", "prop_cases_train",
             "ES2", "MCSE","ES2_thin", "MCSE_thin")
filename <- "output_LOO_nogamma_rho"
if(smote_bool == TRUE){
  filename <- paste0(filename, "_SMOTE")
}
if(interaction_bool == FALSE){
  filename <- paste0(filename, "_nointeractions")
}
if(neigh_bool == FALSE){
  filename <- paste0(filename, "_noneigh")
}
if(reweight_bool == TRUE){
  filename <- paste0(filename,ifelse(alpha_vs_beta_bool, "_a-reweight", "_b-reweight"))
}

if(noRE_bool == TRUE){
  filename <- paste0(filename, "_noRE")
}
filename <- paste0(filename, "_newcov_", outcome, "_",prior,".RData")
save(list = varlist, file = paste0(fldr, "scripts/", filename))

