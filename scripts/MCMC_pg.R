require(MASS)
require(mvtnorm)
require(BayesLogit)
require(mvtnorm)
inverse_link <- function(alpha){
  1/(1+exp(-(alpha)))
}
get_prior_prec <- function(w.sym, rho, nalpha){
  rho*(diag(rowSums(w.sym))-w.sym) + (1-rho)*diag(nalpha)
}

mcmc <- function(y, mapping, X, niter, w.sym, rho = 0.9, a_rho = 10, b_rho = 10,
                 tau2_beta = 10, tau2_alpha = 1, sample_rho = TRUE, sample_tau = TRUE){
  ## mapping is a dataframe with a column "tract" that tells us which tract (from 1 to n) the observation belongs to
  p1 <- ncol(X)
  alpha_un <- unique(mapping$tract)
  nalpha <- length(alpha_un)
  
  shape_a <- 0.1
  rate_a <- 0.1
  
  shape_b <- 0.1
  rate_b <- 0.1
  
  ## if b0 is a scalar:
  # b0 <- 0          ## b0 ~ N(0, tau2_b0)
  ## if b0 is a vector:
  b0 <- rep(0, p1)   ## b0 ~ N(0, tau2_b0 * I)
  tau2_b0 <- 10^2
  
  a0 <- 0          ## a0 ~ N(0, tau2_a0)
  tau2_a0 <- 10^2
  
  # beta is  N(b0, tau2_beta * I)
  # alpha is N(a0, tau2_alpha * prior_prec_alpha^-1)
  
  if(all(X[,1]==1)){
    warning("mcmc: including a vector of 1's in X introduces identifiability problems with random effects.")
  }
  
  # given mapping, we create the matrix
  mappingMAT <- matrix(0, nrow = nalpha, ncol = length(y))
  for(i in 1:nalpha){
    index <- which(mapping$tract == i)
    mappingMAT[i, index] <- 1
  }
  
  ALPHA <- matrix(0, nalpha, 1+niter)
  BETA <- matrix(0, p1, 1+niter)
  TAU2_ALPHA <- numeric(1+niter)
  TAU2_BETA <- numeric(1+niter)
  ## if b0 is a scalar:
  # B0 <- numeric(1+niter)
  ## if b0 is a vector:
  B0 <- matrix(0, p1, 1+niter)
  A0 <- numeric(1+niter)
  RHO <- numeric(1+niter)
  
  # starting points:
  alpha_cur <- matrix(0, nalpha, 1)
  beta_cur <- matrix(0, p1, 1)
  ALPHA[,1] <- as.numeric(alpha_cur) ## neighborhood fixed effects
  BETA[,1] <- as.numeric(beta_cur)   ## effects of time varying covariates
  TAU2_ALPHA[1] <- tau2_alpha
  TAU2_BETA[1] <- tau2_beta
  ## if b0 is a scalar:
  # B0[1] <- b0
  ## if b0 is a vector:
  B0[,1] <- as.numeric(b0)
  A0[1] <- a0
  RHO[1] <- rho
  
  nu_tau2alpha = nalpha - 1
  
  for(k_iter in 1:niter){
    if(k_iter %% 10 == 0){
      cat(k_iter, " ")
    }
    # cat(k_iter, " ")
    prior_prec_alpha = get_prior_prec(w.sym, rho, nalpha)
    ## sample PG (polya gamma)
    ## tmp_par should be logit(p_i), ie beta * x_i, in our case we also add the random effects
    tmp_par <- as.numeric(alpha_cur[mapping$tract] + X %*% beta_cur) 
    w <- rpg(length(tmp_par), 1, tmp_par)
    y_tilde <- (y - 0.5)/w
    
    if(nalpha == 1){
      prec_post_alpha <- as.numeric(mappingMAT %*% w) + prior_prec_alpha / tau2_alpha
    } else {
      prec_post_alpha <- diag(as.numeric(mappingMAT %*% w)) + prior_prec_alpha / tau2_alpha
    }
    cov_post_alpha <- solve(prec_post_alpha)
    m_alpha_tmp <- rowSums(prior_prec_alpha) * a0 /tau2_alpha + 
      mappingMAT %*% (w * as.numeric(y_tilde - X %*% beta_cur))
    m_alpha <- cov_post_alpha %*% m_alpha_tmp
    alpha_cur <- t(rmvnorm(1, m_alpha, cov_post_alpha))
    
    XtOmega <- t(X * w)
    V_beta <- solve((XtOmega %*% X + (diag(p1) / tau2_beta))) # XtOmega %*% X = t(X) %*% Omega %*% X
    
    m_beta_tmp <- XtOmega %*% (y_tilde - alpha_cur[mapping$tract]) + b0 / tau2_beta
    m_beta <- V_beta %*% m_beta_tmp 
    beta_cur <- t(rmvnorm(1, m_beta, V_beta))
    
    ## if b0 is a scalar:
    # b0 <- rnorm(1, (p1/tau2_beta*mean(beta_cur))/(p1/tau2_beta + 1/tau2_b0), sd = 1/sqrt(p1/tau2_beta + 1/tau2_b0))
    ## if b0 is a vector:
    b0 <- rnorm(p1, (beta_cur/tau2_beta)/(1/tau2_beta + 1/tau2_b0), sd = 1/sqrt(1/tau2_beta + 1/tau2_b0))
    a0 <- rnorm(1, (nalpha/tau2_alpha*mean(alpha_cur))/(nalpha/tau2_alpha + 1/tau2_a0), sd = 1/sqrt(nalpha/tau2_alpha + 1/tau2_a0))
    
    if(sample_tau){
      ## sample tau2_beta
      tau2_beta <- 1/rgamma(1, shape = shape_b + p1/2, rate = rate_b + sum((beta_cur - b0)^2)/2)
      ## sample tau2_alpha - conditional on the value of alpha_cur, but not on gamma_cur
      res <- alpha_cur - a0
      s2_tau2alpha <- t(res) %*% prior_prec_alpha %*% res
      tau2_alpha <- 1/rgamma(1, shape = shape_a + nu_tau2alpha/2, rate = rate_a + s2_tau2alpha/2)
    }
    
    if(sample_rho){
      # my proposal is beta(a,b) with a = b*mu/(1-mu) and b = 10
      b <- 5
      rho_new <- rbeta(1, shape1 = b*rho/(1-rho), shape2 = b)
      
      prior_prec_new <- get_prior_prec(w.sym, rho_new, nalpha)
      
      res <- alpha_cur - a0
      post_ratio <- as.numeric(-t(res)%*%(prior_prec_new-prior_prec_alpha)%*%(res)/(2*tau2_alpha))
      post_ratio <- post_ratio + sum(log(diag(chol(prior_prec_new)))-log(diag(chol(prior_prec_alpha))))
      
      prior_ratio <- log(dbeta(rho_new, a_rho,b_rho))-log(dbeta(rho, a_rho,b_rho))
      prop_ratio <- log(dbeta(rho, shape1 = b*rho_new/(1-rho_new), shape2 = b)) - log(dbeta(rho_new, shape1 = b*rho/(1-rho), shape2 = b))
      
      logp_acc <- min(0, prop_ratio+prior_ratio+post_ratio)
      if (log(runif(1)) <= logp_acc) {
        rho <- rho_new
      }
      
    }
    
    ALPHA[, 1+k_iter] <- as.numeric(alpha_cur)
    BETA[, 1+k_iter] <- as.numeric(beta_cur)
    TAU2_ALPHA[1+k_iter] <- tau2_alpha
    TAU2_BETA[1+k_iter] <- tau2_beta
    ## if b0 is a scalar:
    # B0[1+k_iter] <- b0
    ## if b0 is a vector:
    B0[,1+k_iter] <- as.numeric(b0)
    A0[1+k_iter] <- a0
    RHO[1+k_iter] <- rho
  }
  
  return(list(alpha = ALPHA, beta = BETA, rho = RHO,
              tau2_alpha = TAU2_ALPHA, tau2_beta = TAU2_BETA, 
              b0 = B0, a0 = A0))
}

## Check it works on a small dataset
# wtmp <- matrix(c(0,1,1,0,0,
#                  1,0,1,1,0,
#                  1,1,0,0,1,
#                  0,1,0,0,1,
#                  0,0,1,1,0), ncol = 5)
# rho <- 0.9
# nalpha <- nrow(wtmp)
# prior_prec = rho*(diag(rowSums(wtmp))-wtmp) + (1-rho)*diag(nalpha)
# n_tr <- 5
# S <- 1000
# p1 <- 5
# mapping <- data.frame(tract = rep(1:n_tr, each = S),row_index= 1:(S*n_tr))
# alpha <- rmvnorm(1, sigma = prior_prec*2)
# # alpha <- c(-2,-1,0,1,2)
# X <- matrix(rnorm(n_tr*S*p1), ncol = p1)
# beta <- c(1,2,-1, -4, 10)
# theta <- inverse_link( X %*% beta + alpha[mapping$tract] )
# Y <- rbinom(n_tr*S, size = 1, prob = theta)
# mcmc_niter <- NITER <- 1000
# tmp <- mcmc(y = Y,
#             mapping = mapping,
#             X = X,
#             niter = NITER,
#             w.sym = wtmp,
#             tau2_alpha = 0.1, tau2_beta = 0.1)