##################################################################
################### create bayesp_illustration ###################
##################################################################

png("bayesp_illustration.png",width = 10, height = 5, units = "in", res = 300)
par(mfrow =c(1,2),oma = c(0, 0, 2, 0))
x <- seq(-3.5,3.5, length = 100)
p1 <- dnorm(x, mean = 1, sd = 0.5)
plot(x, p1, type = "l", lwd = 2, ylab = "p(x)", ylim = c(0,0.8))
abline(h = 0)
polygon(c(0,x[x>=0]), c(0,p1[x>=0]), col = "gray")
axis(1)
text(0.5, 0.7, labels = paste0("Bayes-p = ",round(1-pnorm(0, mean = 1,sd = 0.5),3)), cex = 1.2, adj = 1)

p2 <- dnorm(x, mean = -.5, sd = 0.8)
plot(x, p2, type = "l", lwd = 2, ylab = "p(x)", ylim = c(0,0.8))
abline(h = 0)
polygon(c(0,x[x<=0]), c(0,p2[x<=0]), col = "gray")
axis(1)
text(0, 0.5, labels = paste0("Bayes-p = ",round(pnorm(0, mean = -.5, sd = 0.8),3)), cex = 1.2, adj = 0)
mtext("Bayes-p illustration", outer = TRUE, cex = 1.5, line = -2)
dev.off()

##################################################################
##########################  rho plots   ##########################
##################################################################

burnin <- 500
thin <- 5
mcmc_niter <- 1000*thin+burnin 
index_thinning <- seq(burnin, mcmc_niter, by = thin)

year <- 8; wdstr <- "results/"; add_str2 <- "_YEAR8"
add_str = "" 
for(output_string in c("PRETERM", "STILLBIRTH")){
  tmp <-load(paste0(wdstr,"output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
  year <- 8; output1 <- output_list1[[year]]; output2 <- output_list2[[year]]
  rhos <- c(output1$rho[index_thinning],output2$rho[index_thinning])
  assign(paste0("rhos_",output_string),rhos)
}

library(bde)

png(filename = "rho_density.png", width = 10, height = 5, units = "in", res = 300)
par(mfrow= c(1,2))
plot(bde(rhos_STILLBIRTH,estimator="boundarykernel",
         lower.limit = 0,upper.limit = 1,options = list(mu=1)), 
     ylab = "", main = "Stillbirth", xlab = expression(rho))
plot(bde(rhos_PRETERM,estimator="boundarykernel",
         lower.limit = 0,upper.limit = 1,options = list(mu=1)), 
     ylab = "", main = "Preterm Birth", xlab = expression(rho))
dev.off()

##################################################################
####################### plot logodds (CAR) #######################
##################################################################

wdstr <- "results/"

bayesp <- function(x){max(mean(x>0),mean(x<0))}
var_neighborhood <- c("proportion Asian","proportion Hispanic","proportion Black", # "Prop_White", 
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
thin <- 5
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

year <- 8; wdstr <- "results/"; add_str2 <- "_YEAR8"
for(output_string in c("PRETERM", "STILLBIRTH")){
  add_str = ""
  
  tmp <-load(paste0(wdstr,"output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
  year <- 8; output1 <- output_list1[[year]]; output2 <- output_list2[[year]]
  beta_tr <- cbind(output1$beta[,index_thinning, drop = FALSE],
                   output2$beta[,index_thinning, drop = FALSE])
  scale_sds <- scale_sds_list[[year]]
  for(i in col_to_be_scaled){
    beta_tr[i,] <- beta_tr[i,] / scale_sds[i]
  }
  beta_LOR <- rbind(rowMeans(beta_tr[,index]),
                    apply(beta_tr[,index], MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975)))
  beta_LOR <- t(beta_LOR[c(2,1,3),])
  colnames(beta_LOR)[2] <- "Mean"; rownames(beta_LOR) <- c(var_individual, var_neighborhood)
  beta_LOR <- beta_LOR[rev(index),]
  
  res <- beta_LOR
  assign(ifelse(output_string == "PRETERM", "res_pre","res_still"),res)
}


### Let's make the plot
filename_str <- "cov_coeff_comb.png"

## let's add an empty line
res <- res_still
res2 <- res[1:14,]
res2 = rbind(res2,NA)
res2 = rbind(res2,res[15:19,])
res_st = res2

res <- res_pre
res2 <- res[1:14,]
res2 = rbind(res2,NA)
res2 = rbind(res2,res[15:19,])
res_pr = res2


k <- nrow(res_st)
labels <- rownames(res_st)

# parameters for plot
las <- 1
length <- 0
angle <- 30
code <- 3
pchs <- rep(19, length.out = k)
space_above <- space_below <- 0.5


png(filename = filename_str, width = 10, height = 7, units = "in", res = 100)

laymat <- matrix(c(1,2,1,3), ncol = 2, nrow = 2)
layout(mat = laymat, heights = c(0.1,2), widths = c(2,1))
# layout.show(3)

par(mar=c(0,0,1,0),mgp=c(0,0,0))
plot(0, 0, xlab = "", ylab = "",
     axes = FALSE, type = "n")  
mtext("Log odds-ratio", side =3, cex = 1.5, line = -1.2, font = 2,
      at=par("usr")[1]+0.65*diff(par("usr")[1:2]))

par(mar=c(3,14,2,-0.1)+0.1,mgp=c(5,1,0), cex.axis = 1.3)

res <- res_st
xlim <- range(res, na.rm = T); xlim[1] <- xlim[1] - 0.3
ylim <- c(1 - space_below, k + space_above)

plot(0, 0, xlim = xlim, ylim = ylim, xlab = "", ylab = "",
     axes = FALSE, type = "n", las = las)    
abline(v = 0, lty = 2, lwd = 1)
axis(1)
axis(2, at = (1:k)[-15], labels = labels[-15], las = las)
box()
title("Stillbirth", cex = 0.8, line = 0.5)
offset = +0.15
col <- rep(1, length.out = k)
arrows(res[,c(1)], 1:k + offset,res[,c(3)], 1:k + offset, lty = 1, lwd = 1, col = col,
       length = length, angle = angle, code = code)
points(res[,c(2)], 1:k + offset, pch = pchs, col = col)

par(mar=c(3,-0.1,2,1)+0.1,mgp=c(5,1,0))
res <- res_pr
xlim <- range(res, na.rm = T); xlim[1] <- xlim[1] - 0.3
ylim <- c(1 - space_below, k + space_above)

plot(0, 0, xlim = xlim, ylim = ylim, xlab = "", ylab = "",
     axes = FALSE, type = "n", las = las)
abline(v = 0, lty = 2, lwd = 1)
axis(1)
box()
title("Preterm birth", cex = 0.8, line = 0.5)
offset = +0.15
col <- rep(1, length.out = k)
arrows(res[,c(1)], 1:k + offset,res[,c(3)], 1:k + offset, lty = 1, lwd = 1, col = col,
       length = length, angle = angle, code = code)
points(res[,c(2)], 1:k + offset, pch = pchs, col = col)

dev.off()

##################################################################
################### Neighborhood cluster risk ####################
##################################################################

## See neigh_cluster_analysis.R for this plot






### The following are for the supplement

#############################################################################
######################## MCMC diagnostics (figure S4) #######################
#############################################################################

var_neighborhood <- c("proportion Asian","proportion Hispanic","proportion Black", # "Prop_White", 
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
  
  tmp1 <- output_list1[[year]]; tmp2 <- output_list2[[year]]
  
  assign(paste0("tmp1_",output_string),tmp1)
  assign(paste0("tmp2_",output_string),tmp2)
}

# STILLBIRTH

output1 <- tmp1_STILLBIRTH
output2 <- tmp2_STILLBIRTH

betas <- cbind(output1$beta[,index_thinning, drop = FALSE],
               output2$beta[,index_thinning, drop = FALSE])
png(filename = "traceplots_betas_noSMOTE_STI_car.png", width = 10, height = 5, units = "in", res = 300)
is <- c(4, 11, 15)

par(mfrow = c(2,3))
i <- is[1]
plot(output1$beta[i,], type = "l", col = rgb(0,0,0,0.5), lwd = 1, ylab = "", main = paste("Trace plot:",vars[i]))
lines(output2$beta[i,], col = rgb(1,0,0,0.5), lwd = 1)

i <- is[2]
plot(output1$beta[i,], type = "l", col = rgb(0,0,0,0.5), lwd = 1,  ylab = "", main = paste("Trace plot:",vars[i]))
lines(output2$beta[i,], col = rgb(1,0,0,0.5), lwd = 1)

i <- is[3]
plot(output1$beta[i,], type = "l", col = rgb(0,0,0,0.5), lwd = 1,  ylab = "", main = paste("Trace plot:",vars[i]))
lines(output2$beta[i,], col = rgb(1,0,0,0.5), lwd = 1)

acf(betas[is[1],], ylab = "", main = paste("ACF:",vars[is[1]]))
acf(betas[is[2],], ylab = "", main = paste("ACF:",vars[is[2]]))
acf(betas[is[3],], ylab = "", main = paste("ACF:",vars[is[3]]))
dev.off()

# PRETERM

output1 <- tmp1_PRETERM
output2 <- tmp2_PRETERM
png(filename = "traceplots_betas_noSMOTE_PRE_ind.png", width = 10, height = 5, units = "in", res = 300)
is <- c(1, 6, 18)

par(mfrow = c(2,3))
i <- is[1]
plot(output1$beta[i,], type = "l", col = rgb(0,0,0,0.5), lwd = 1, ylab = "", main = paste("Trace plot:",vars[i]))
lines(output2$beta[i,], col = rgb(1,0,0,0.5), lwd = 1)

i <- is[2]
plot(output1$beta[i,], type = "l", col = rgb(0,0,0,0.5), lwd = 1,  ylab = "", main = paste("Trace plot:",vars[i]))
lines(output2$beta[i,], col = rgb(1,0,0,0.5), lwd = 1)

i <- is[3]
plot(output1$beta[i,], type = "l", col = rgb(0,0,0,0.5), lwd = 1,  ylab = "", main = paste("Trace plot:",vars[i]))
lines(output2$beta[i,], col = rgb(1,0,0,0.5), lwd = 1)

acf(betas[is[1],], ylab = "", main = paste("ACF:",vars[is[1]]))
acf(betas[is[2],], ylab = "", main = paste("ACF:",vars[is[2]]))
acf(betas[is[3],], ylab = "", main = paste("ACF:",vars[is[3]]))
dev.off()


