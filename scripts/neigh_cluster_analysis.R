library(sp)
library(rgeos) # needed for maptools
library(ggplot2)
library(ggmap)
library(dplyr)
library(maptools)
library(xtable)

source("scripts/functions_scale.R")
load("data/phillytracts")
load("data/phillygooglemap.rdata")
wdstr <- "results/"


#######################################
########## Find the clusters ##########
#######################################

######## Extract the data first #######


burnin <- 500
thin <- 5
mcmc_niter <- 1000*thin+burnin 
index_thinning <- seq(burnin, mcmc_niter, by = thin)

var_neighborhood <- c("prAsian", "prHispanic", "prBlack", # "Prop_White", 
                      "women", "poverty", 
                      "public assistance", "labor force", 
                      "birth", "high school", 
                      "college", "occupied housing", "violation", 
                      "violent", "nonviolent")
var_individual <- c("Hispanic","Black", "Asian", "multiple birth", "age")
p1 <- length(var_individual)
p2 <- length(var_neighborhood)
col_to_be_scaled <- p1:(p1+p2)

var_str <- c("avg_phat_samples", "avg_yhat_samples","avg_phat2_samples", 
             "avg_phat2shift_samples","avg_phat_shift_samples")

year <- 8; add_str = ""; add_str2 <- "_YEAR8"
for(output_string in c("PRETERM", "STILLBIRTH")){
  tmp <-load(paste0(wdstr,"output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
  scale_sds <- scale_sds_list[[year]]
  scale_mus <- scale_means_list[[year]]
  
  output1 <- output_list1[[year]]
  output2 <- output_list2[[year]]
  alphas <- cbind(output1$alpha[,index_thinning, drop = FALSE],
                  output2$alpha[,index_thinning, drop = FALSE])
  alphas <- rowMeans(alphas)
  betas <- cbind(output1$beta[,index_thinning, drop = FALSE],
                 output2$beta[,index_thinning, drop = FALSE])
  beta <- rowMeans(betas)
  beta <- beta[col_to_be_scaled]/as.numeric(scale_sds[col_to_be_scaled])
  alphas <- alphas - sum(beta * as.numeric(scale_mus[col_to_be_scaled]))
  
  data_CAR <- data.frame(tractID = tractID_unique_list[[year]],  
                         phat_all = mean_phat_all_list[[year]], phat2_all = mean_phat2_all_list[[year]],
                         phat2shift_all = mean_phat2shift_all_list[[year]], phat_neig = phat_neig_list[[year]],
                         phatshift_neigh = phat_neigshift_list[[year]],
                         alpha = alphas)

  for(x in var_str){
    assign(ifelse(output_string == "PRETERM",paste0(x,"_PRE"),paste0(x,"_STILL")), get(x))
  }
  for(x in var_str){
    data_CAR <- cbind(data_CAR, rowMeans(get(x)))
    names(data_CAR)[ncol(data_CAR)] <- x
  }

  assign(ifelse(output_string == "PRETERM","data_CAR_PRE","data_CAR_STILL"),data_CAR)
}

phat_bool <- TRUE     # we are using phat instead of phat2
shift_bool <- FALSE   # shift (i.e. calibration) is needed only for reweight and SMOTE

########## Find the clusters ##########

var <- ifelse(phat_bool, ifelse(shift_bool, "avg_phat_shift_samples", "avg_phat_samples"), 
              ifelse(shift_bool, "avg_phat2shift_samples", "avg_phat2_samples"))
var_orig_PRE <- data_CAR_PRE[,var]
var_orig_STILL <- data_CAR_STILL[,var]

var_samples_str <- paste0(ifelse(phat_bool,ifelse(shift_bool,"avg_phat_shift_samples","avg_phat_samples"),
                                 ifelse(shift_bool,"avg_phat2shift_samples","avg_phat2_samples")))
var_samples_PRE <- get(paste0(var_samples_str, "_PRE"))
var_samples_STILL <- get(paste0(var_samples_str, "_STILL"))

set.seed(123456789)
km_PRE <- kmeans(var_orig_PRE, centers = 3, nstart = 100)
km_STILL <- kmeans(var_orig_STILL, centers = 3, nstart = 100)
var1 <- km_PRE$centers[km_PRE$cluster]
var2 <- km_STILL$centers[km_STILL$cluster]

tmp_order <- order(km_PRE$centers)
cl1 <- match(km_PRE$cluster, tmp_order)

tmp_order <- order(km_STILL$centers)
cl2 <- match(km_STILL$cluster, tmp_order)

#################################################
######## Density of cluster prob (Fig 6) ########
#################################################

for(outcome in c("PRETERM","STILLBIRTH")){
  var_samples <- get(ifelse(outcome == "PRETERM","var_samples_PRE","var_samples_STILL"))
  for(x in 1:3){
    ps <- colMeans(var_samples[cl1 == x,])
    assign(paste0("ps",x),ps)
  }
  
  filename_str <- paste0("density_cl_phat", ifelse(phat_bool,"","2"),ifelse(shift_bool,"_shift",""))
  filename_str <- paste0(filename_str,ifelse(outcome == "PRETERM", "_PRE","_STILL"),".png")
  
  xlims <- range(c(density(ps1)$x, density(ps2)$x,density(ps3)$x))
  ylims <- range(c(density(ps1)$y, density(ps2)$y,density(ps3)$y))
  
  png(filename = filename_str, width = 6, height = 3.5, units = "in", res = 100)
  par(mar = c(4.1,2.6,4.1,2.1))
  plot(density(ps1), xlim = xlims, ylim = ylims * 1.2, ylab = "",
       xlab = "Risk probability", lty = 2, lwd = 2, col = "springgreen4",
       main = paste0(ifelse(outcome == "PRETERM","Preterm birth ","Stillbirth "),"neighborhood cluster risks"))
  lines(density(ps2), lwd = 2, lty = 1, col = "orange")
  lines(density(ps3), lwd = 2, lty = 3, col = "red3")
  legend("topright", c("Lower risk", "Moderate risk", "Higher risk"),
         lty = c(2,1,3), lwd = rep(2,3), col = c("springgreen4", "orange", "red3"), ncol = 3, bty ="n")
  dev.off()
  
  CI <- rbind(quantile(ps1, probs = c(0.025,0.975)),
              quantile(ps2, probs = c(0.025,0.975)),
              quantile(ps3, probs = c(0.025,0.975)))
  assign(paste0("CI_",ifelse(outcome == "PRETERM","PRE","STILL")),CI)
}

CI_PRE
CI_STILL

#################################################
##########  Create CSV for supplement ###########
#################################################

year <- 8; add_str = ""; add_str2 <- "_YEAR8"
for(output_string in c("PRETERM", "STILLBIRTH")){
  tmp <-load(paste0(wdstr,"output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
  
  data_CAR <- get(paste0("data_CAR_",ifelse(output_string == "PRETERM","PRE","STILL")))
  cl <- get(paste0("cl", ifelse(output_string == "PRETERM", "1","2")))
  data_CSV <- data_CAR[,c("tractID","avg_phat_samples")]
  data_CSV <- cbind( data_CSV, 
                t(apply(avg_phat_samples, MARGIN = 1, quantile, probs = c(0.025,0.975))),
                cl)
  
  assign(ifelse(output_string == "PRETERM","data_CSV_PRE","data_CSV_STILL"),data_CSV)
}
data_CSV <- cbind(data_CSV_PRE, data_CSV_STILL)
data_CSV <- data_CSV[,c(1,10,5,7,2,8,9,3,4)]
colnames(data_CSV) <- c("tractID", "risk_cat_stillbirth","risk_cat_preterm",
                        "risk_prob_stillbirth","risk_prob_preterm",
                        "lCI_prob_stillbirth","uCI_prob_stillbirth",
                        "lCI_prob_preterm","uCI_prob_preterm")

tmpf <- function(x){
  ifelse(x == 1, "Lower", ifelse(x == 2, "Moderate", "Higher"))
}

data_CSV[,2] <- tmpf(data_CSV[,2])
data_CSV[,3] <- tmpf(data_CSV[,3])
data_CSV[,4:9] <- round(data_CSV[,4:9], 4)
head(data_CSV)

write.table(data_CSV, file = "Neighborhood_probabilities.csv", sep = ",", row.names = F)

#################################################
### Map of clusters (bottom panels of Fig 6) ####
#################################################

# Run twice (one for each outcome)
PRE_bool <-  TRUE # FALSE # TRUE

if(PRE_bool){
  var <- var1; group_ind <- match(data_CAR_PRE$tractID, tracts@data$GEOID10)
  vars <- c(var,var_orig_PRE) 
} else {
  var <- var2; group_ind <- match(data_CAR_STILL$tractID, tracts@data$GEOID10)
  vars <- c(var,var_orig_STILL) 
}

polyfortified <- fortify(tracts)
limits <- c(min(vars, na.rm=T), max(vars, na.rm=T))
MI <- limits[1]
MA <- limits[2]
mi <- min(var, na.rm=T)
ma <- max(var, na.rm=T)
z <- (mean(var, na.rm=T) -mi)/(ma-mi)
z_adj<-(ma-mi)/(MA-MI) * z + (mi-MI)/(MA-MI)

my.data <- data.frame(id = as.numeric(rownames(tracts@data[group_ind,])), value =var)
my.data$id <- as.character(my.data$id)
plotData <- left_join(polyfortified,my.data, by = "id")

title_string <- paste0(ifelse(PRE_bool, "Preterm birth", "Stillbirth"), ": neighborhood clusters")
palette <- ifelse(PRE_bool, "RdYlGn", "RdYlBu")

p <- ggmap(googlemap) #ggplot()
p <- p + geom_polygon(data=plotData, aes(x=long,y=lat,group=group,fill=value),
                      color=alpha("black",0.15), alpha =1) + ggtitle(title_string)
p <- p + scale_fill_distiller(palette = palette, limits = c(MI,MA), values = c(0,z_adj,1), name =expression(hat(p)[i]))
p <- p + theme(plot.title = element_text(hjust = 0.5, vjust = -0.5,size=23),
               panel.border =  element_rect(colour = 'transparent', fill = 'transparent'),
               axis.ticks = element_blank(),
               panel.spacing = unit(0, 'mm'), axis.text = element_blank(),
               panel.grid = element_blank(), axis.title = element_blank(),
               legend.key.height = unit(1,"cm"),
               legend.text=element_text(size=15),
               legend.title = element_text(size=15),
               legend.justification = c(1,0.035), legend.position = c(1,0.035),
               legend.background = element_rect(fill= alpha("white",0.6)))
p <- p + scale_bar(lon = -75.17, lat = 39.825, 
                   distance_lon = 5, distance_lat = 0.5, distance_legend = 1.1, 
                   dist_unit = "km", orientation = FALSE)
p

file_name <- paste0("new_", ifelse(phat_bool, ifelse(shift_bool, "phat_shift", "phat"), 
                                   ifelse(shift_bool, "phat2shift", "phat2")),"_km")
file_name <- paste0(file_name, "_", ifelse(PRE_bool, "PRE", "STILL"),"_")
file_name <- paste0(file_name, "noSMOTE_")
file_name <- paste0(file_name, "CAR.png")
ggsave(file_name, p, scale = 1, width = 8, height = 8, units = 'in')


#############################################################
##### Neighborhood cluster analysis (Tables 3, S2, S3) ######
#############################################################

################### Extract the data first ##################

philly_all_covariates_med <- as.data.frame(philly_all_covariates %>% group_by(tractID) %>% summarize_all(median))

##### ANOVA

Fvalues_PRE <- numeric(p2); names(Fvalues_PRE) <- var_neighborhood_old
pvalues_PRE <- numeric(p2); names(pvalues_PRE) <- var_neighborhood_old
Fvalues_STILL <- numeric(p2); names(Fvalues_STILL) <- var_neighborhood_old
pvalues_STILL <- numeric(p2); names(pvalues_STILL) <- var_neighborhood_old
for(col in var_neighborhood_old){
  m<- aov(philly_all_covariates_med[,col] ~ as.factor(km_PRE$cluster))
  s <- summary(m)
  Fvalues_PRE[col] <- s[[1]][1,4]
  pvalues_PRE[col] <- s[[1]][1,5]
  m<- aov(philly_all_covariates_med[,col] ~ as.factor(km_STILL$cluster))
  s <- summary(m)
  Fvalues_STILL[col] <- s[[1]][1,4]
  pvalues_STILL[col] <- s[[1]][1,5]
}
names(Fvalues_PRE) <- var_neighborhood; names(pvalues_PRE) <- var_neighborhood
names(Fvalues_STILL) <- var_neighborhood; names(pvalues_STILL) <- var_neighborhood

##### Regression

cl_phat_regr_bool <- TRUE

beta_PRE <- numeric(p2); names(beta_PRE) <- var_neighborhood_old
b_pvalues_PRE <- numeric(p2); names(b_pvalues_PRE) <- var_neighborhood_old
beta_STILL <- numeric(p2); names(beta_STILL) <- var_neighborhood_old
b_pvalues_STILL <- numeric(p2); names(b_pvalues_STILL) <- var_neighborhood_old
for(col in var_neighborhood_old){
  if(cl_phat_regr_bool){
    m <- lm(philly_all_covariates_med[,col] ~ var1)
  } else {
    m <- lm(philly_all_covariates_med[,col] ~ cl1)
  }
  s <- summary(m)
  beta_PRE[col] <- s$coefficients[2,1]
  b_pvalues_PRE[col] <- s$coefficients[2,4]
  if(cl_phat_regr_bool){
    m <- lm(philly_all_covariates_med[,col] ~ var2)
  } else {
    m <- lm(philly_all_covariates_med[,col] ~ cl2)
  }
  s <- summary(m)
  beta_STILL[col] <- s$coefficients[2,1]
  b_pvalues_STILL[col] <- s$coefficients[2,4]
}
names(beta_PRE) <- var_neighborhood; names(b_pvalues_PRE) <- var_neighborhood
names(beta_STILL) <- var_neighborhood; names(b_pvalues_STILL) <- var_neighborhood

################## Create table 3 #################

for(PRE_bool in c(TRUE, FALSE)){
	if(PRE_bool){
	  phat_neig <- var_orig_PRE # using the same data we had to create clusters
	  km_tmp <- km_PRE
	  str_append <- "_PRE"
	  cl_tmp <- cl1
	} else {
	  phat_neig <- var_orig_STILL # using the same data we had to create clusters
	  km_tmp <- km_STILL
	  str_append <- "_STILL"
	  cl_tmp <- cl2
	}

	i <- 1
	index_i <- which(cl_tmp == i)

	dataset <- data.frame()
	dataset <- rbind(dataset, as.numeric(apply(philly_all_covariates_med[index_i,],2,mean)))
	colnames(dataset) <- colnames(philly_all_covariates)
	dataset$phat <- mean(phat_neig[index_i])
	dataset$count <- length(index_i)
	for(i in 2:3){
	  index_i <- which(cl_tmp == i)
	  tmp <- as.numeric(apply(philly_all_covariates_med[index_i,],2,mean))
	  tmp <- c(tmp, mean(phat_neig[index_i]), length(index_i))
	  names(tmp) <- c(colnames(philly_all_covariates),"phat")
	  dataset <- rbind(dataset, tmp)
	}

	data_print <- dataset[,-c(1,2)]
	colnames(data_print)[match(var_neighborhood_old,colnames(data_print))] <- var_neighborhood
	data_print[, var_neighborhood] <- round(data_print[, var_neighborhood],2) 
	data_print$phat <- round(100*data_print$phat,2)
	data_print$phat <- paste0(data_print$phat,"%")
	data_print <- data_print[, c(15,16,1:14)]

	#### Include only + / - / o ####
	plusminus <- ifelse(get(paste0("beta", str_append))>0,"+","--")
	plusminus <- sapply(1:p2, FUN = function(i) ifelse(get(paste0("b_pvalues", str_append))[i] < 0.05, plusminus[i], " ")) 

	data_print2 <- t(data_print)
	data_print2 <- cbind(data_print2,c(NA,NA, plusminus))
	assign(paste0("data_print2",str_append),data_print2)
	assign(paste0("dataset",str_append),dataset)
}

data_print2_comb <- cbind(data_print2_STILL,data_print2_PRE)
xtable(data_print2_comb) ### Table 3 ###

################## Create table S2,S3 #################

### PRETERM
PRE_bool = TRUE
dataset <- get(paste0("dataset",str_append))

data_print <- dataset[,-c(1,2)]
colnames(data_print)[match(var_neighborhood_old,colnames(data_print))] <- var_neighborhood
data_print$phat <- round(100*data_print$phat,2)
data_print$phat <- paste0(data_print$phat,"%")
data_print <- data_print[, c(15,1:14)] # we can remove phat and count

data_print3 <- rbind(data_print,c(NA, get(paste0("Fvalues", str_append))),
                                c(NA, get(paste0("pvalues", str_append))))
data_print3 <- rbind(data_print3, 
                     c(NA, get(paste0("beta", str_append))), 
                     c(NA, get(paste0("b_pvalues", str_append))))
data_print3[,-1] <- round(data_print3[,-1],4)
data_print3 <- t(data_print3)
data_print3 <- data_print3[-1,] # we can actually remove phat!
xtable(data_print3)


### STILLBIRTH
PRE_bool = FALSE
dataset <- get(paste0("dataset",str_append))

data_print <- dataset[,-c(1,2)]
colnames(data_print)[match(var_neighborhood_old,colnames(data_print))] <- var_neighborhood
data_print$phat <- round(100*data_print$phat,2)
data_print$phat <- paste0(data_print$phat,"%")
data_print <- data_print[, c(15,1:14)] # we can remove phat and count

data_print3 <- rbind(data_print,c(NA, get(paste0("Fvalues", str_append))),
                                c(NA, get(paste0("pvalues", str_append))))
data_print3 <- rbind(data_print3, 
                     c(NA, get(paste0("beta", str_append))), 
                     c(NA, get(paste0("b_pvalues", str_append))))
data_print3[,-1] <- round(data_print3[,-1],4)
data_print3 <- t(data_print3)
data_print3 <- data_print3[-1,] # we can actually remove phat!
xtable(data_print3) 

