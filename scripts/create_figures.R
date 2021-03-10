################### create bayesp_illustration ###################

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

################### plot logodds ###################

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

# order_covariates <- c("age","white","black","Hispanic","Asian","multiple_birth",
#                       "Prop_Asian","Prop_Hispanic","Prop_Black","prop_women_15_to_50",
#                       "prop_women_below_poverty","prop_women_public_assistance",
#                       "prop_women_labor_force","prop_birth_last_12_months",
#                       "prop_women_HS_grad","prop_women_college_grad",
#                       "log_occupied_housing","log_housing_violation",
#                       "log_violent_crime","log_nonviolent_crime")
# index <- match(order_covariates,c(var_individual, var_neighborhood))



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


beta_LOR <- rbind(rowMeans(beta_tr[,index]),
                 apply(beta_tr[,index], MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975)))
beta_LOR <- t(beta_LOR[c(2,1,3),])
colnames(beta_LOR)[2] <- "Mean"
rownames(beta_LOR) <- c(var_individual, var_neighborhood)
# beta_LOR <- beta_LOR[index,]

beta_LOR_noSMOTE <- beta_LOR
beta_LOR_SMOTE <- beta_LOR

res <- cbind(beta_LOR_noSMOTE, beta_LOR_SMOTE)
res <- res[order(res[,2]),]

# res_pre <- res
# res_still <- res

### Let's make the plot

k <- nrow(res)
labels <- rownames(res)

# parameters for plot
las <- 1
length <- 0
angle <- 30
code <- 3
pchs <- rep(19, length.out = k)
space_above <- space_below <- 0.5

xlim <- range(res); xlim[1] <- xlim[1] - 0.3
ylim <- c(1 - space_below, k + space_above)

# png(filename = "cov_coeff_pre.png", width = 7, height = 7, units = "in", res = 300)
png(filename = "cov_coeff_still.png", width = 7, height = 7, units = "in", res = 300)

par(mar=c(3,14,4,1)+0.1,mgp=c(5,1,0))
plot(0, 0, xlim = xlim, ylim = ylim, xlab = "", ylab = "",
     axes = FALSE, type = "n", las = las)    
abline(v = 0, lty = 2, lwd = 1)
axis(1)
axis(2, at = 1:k, labels = labels, las = las)
box()
# legend(-2,14, c("noSMOTE","SMOTE"), col = 1:2, lty = rep(1,2), pch = rep(pchs[1],2), box.lty=1, box.lwd=0.5)
# legend(-3.5,14, c("noSMOTE","SMOTE"), col = 1:2, lty = rep(1,2), pch = rep(pchs[1],2), box.lty=1, box.lwd=0.5)
legend("topleft", c("noSMOTE","SMOTE"), col = 1:2, lty = rep(1,2), pch = rep(pchs[1],2))
title("Log odds-ratio")

offset = +0.15
col <- rep(1, length.out = k)
arrows(res[,c(1)], 1:k + offset,res[,c(3)], 1:k + offset, lty = 1, lwd = 1, col = col,
       length = length, angle = angle, code = code)
points(res[,c(2)], 1:k + offset, pch = pchs, col = col)


offset = -0.15
col <- rep(2, length.out = k)
arrows(res[,4], 1:k + offset, res[,6], 1:k + offset, lty = 1, lwd = 1, col = col,
       length = length, angle = angle, code = code)
points(res[,5], 1:k + offset, pch = pchs, col = col)

dev.off()
