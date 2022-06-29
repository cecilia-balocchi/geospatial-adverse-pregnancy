# rm(list = ls())

library(sp)
library(rgeos) # needed for maptools
library(ggplot2)
library(ggmap)
library(dplyr)
library(maptools)

source("scripts/functions_scale.R")
load("data/phillytracts")
load("data/phillygooglemap.rdata")
wdstr <- "results/"

#############################################################################
#################### Proportion of outcomes map (fig 1) #####################
#############################################################################

tmp = load(paste0(wdstr,"plot_ratios_by_tract.RData"))
var1 <- ratio_by_tract_after_preprocessing$preterm_ratio
var2 <- ratio_by_tract_after_preprocessing$stillbirth_ratio
var3 <- deliveries_by_tract_after_preprocessing$deliveries

group_ind <- match(ratio_by_tract_after_preprocessing$tractID, as.character(tracts@data$GEOID10))
all.equal(deliveries_by_tract_after_preprocessing$tractID, ratio_by_tract_after_preprocessing$tractID)

## change these to make the desired plot. 
## Figure 1, left panel was created with PRE_bool = FALSE; ndel_bool = FALSE
## Figure 1, right panel was created with PRE_bool = TRUE; ndel_bool = FALSE
PRE_bool = TRUE
ndel_bool = FALSE
if(ndel_bool){
  var <- var3
} else {
  if(PRE_bool){
    var <- var1
  } else {
    var <- var2
  }
}
vars <- var

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

if(!ndel_bool){
  title_string <- ifelse(PRE_bool, "Preterm birth", "Stillbirth")
  palette <- ifelse(PRE_bool, "RdYlGn", "RdYlBu")
} else {
  title_string <- "Number of deliveries"
  palette <- "BuPu"
}

p <- ggmap(googlemap) # use  ggplot() if you don't have a map
p <- p + geom_polygon(data=plotData, aes(x=long,y=lat,group=group,fill=value),
                      color=alpha("black",0.15), alpha =1) + ggtitle(title_string)
if(!ndel_bool){
  p <- p + scale_fill_distiller(palette = palette, limits = c(MI,MA),values = c(0,z_adj,1), name =expression(prop[i]))#name =expression(hat(p)[i]))
} else {
  p <- p + scale_fill_distiller(palette = palette, limits = c(MI,MA),values = c(1,z_adj,0), name =expression(n[i]))#name =expression(hat(p)[i]))
}
p <- p + theme(plot.title = element_text(hjust = 0.5, vjust = -1,size=23),
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
file_name <- ifelse(ndel_bool, "n_data", "phat_data")
file_name <- ifelse(ndel_bool,paste0(file_name, ".png"),paste0(file_name, "_", ifelse(PRE_bool, "PRE", "STILL"),".png"))
ggsave(file_name, p, scale = 1, width = 8, height = 8, units = 'in')


#############################################################################
########### Map of predicted probabilities (top panels of fig 5) ############
#############################################################################

burnin <- 500
thin <- 5
mcmc_niter <- 1000*thin+burnin
index_thinning <- seq(burnin, mcmc_niter, by = thin)

bayesp <- function(x){max(mean(x>0),mean(x<0))}
var_neighborhood <- c("prop_Asian", "prop_Hispanic", "prop_Black", # "Prop_White", 
                      "prop_women_15_to_50", "prop_women_below_poverty", 
                      "prop_women_public_assistance", "prop_women_labor_force", 
                      "prop_birth_last_12_months", "prop_women_HS_grad", 
                      "prop_women_college_grad", "log_occupied_housing", "log_housing_violation", 
                      "log_violent_crime", "log_nonviolent_crime")
var_individual <- c("Hispanic","Black", "Asian", "multiple_birth", "age")


p1 <- length(var_individual)
p2 <- length(var_neighborhood)
col_to_be_scaled <- p1:(p1+p2)

var_str <- c("avg_phat_samples", "avg_yhat_samples","avg_phat2_samples", 
             "avg_phat2shift_samples","avg_phat_shift_samples")

year <- 8; add_str <- ""; add_str2 <- "_YEAR8"
for(output_string in c("PRETERM", "STILLBIRTH")){
  tmp <-load(paste0(wdstr,"output_LOO_nogamma_rho_nointeractions_",add_str,"newcov_",output_string,"_CAR",add_str2,".RData"))
  
  for(x in var_str){
    assign(ifelse(output_string == "PRETERM",paste0(x,"_PRE"),paste0(x,"_STILL")), get(x))
  }
  
  x = "tractID_unique"
  assign(ifelse(output_string == "PRETERM",paste0(x,"_PRE"),paste0(x,"_STILL")), tractID_unique_list[[year]])
}

# Run the following twice, once for each outcome, by chaning PRE_bool to TRUE and FALSE.
PRE_bool <- TRUE # FALSE # TRUE
shift_bool <- FALSE # shift (i.e. calibrating probabilities), is needed only for SMOTE or reweight

if(PRE_bool){
  if(shift_bool){
    var1 <- rowMeans(avg_phat_shift_samples_PRE)
    var2 <- rowMeans(avg_phat2shift_samples_PRE)
  } else {
    var1 <- rowMeans(avg_phat_samples_PRE)
    var2 <- rowMeans(avg_phat2_samples_PRE)
  }
  group_ind <- match(tractID_unique_PRE, tracts@data$GEOID10)
} else {
  if(shift_bool){
    var1 <- rowMeans(avg_phat_shift_samples_STILL)
    var2 <- rowMeans(avg_phat2shift_samples_STILL)
  } else {
    var1 <- rowMeans(avg_phat_samples_STILL)
    var2 <- rowMeans(avg_phat2_samples_STILL)
  }
  group_ind <- match(tractID_unique_STILL, tracts@data$GEOID10)
}

phat_bool <- TRUE  ## we ended up using phat instead of phat2, as more interpretable
if(phat_bool){
  var <- var1
} else {
  var <- var2 
}
vars <- c(var1, var2)


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

title_string <- paste0(ifelse(PRE_bool, "Preterm birth", "Stillbirth"),": neighborhood probabilities")
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

file_name <- paste0("new_", ifelse(phat_bool, ifelse(shift_bool, "phat_shift", "phat"), 
                    ifelse(shift_bool, "phat2shift", "phat2")))
file_name <- paste0(file_name, "_", ifelse(PRE_bool, "PRE", "STILL"),"_")
file_name <- paste0(file_name, ifelse(reweight_bool, "reweight", "noSMOTE"),"_")
file_name <- paste0(file_name, ifelse(PRE_bool, "CAR.png", "CAR.png"))
ggsave(file_name, p, scale = 1, width = 8, height = 8, units = 'in')


#############################################################################
########### Map of cluster assignments (bottom panels of fig 5) #############
#############################################################################

## See neigh_cluster_analysis.R for this map


#############################################################################
###################### Map of random effects (Fig S1) #######################
#############################################################################

burnin <- 500
thin <- 5
mcmc_niter <- 1000*thin+burnin
index_thinning <- seq(burnin, mcmc_niter, by = thin)

bayesp <- function(x){max(mean(x>0),mean(x<0))}
var_neighborhood <- c("prop_Asian", "prop_Hispanic", "prop_Black", # "Prop_White", 
                      "prop_women_15_to_50", "prop_women_below_poverty", 
                      "prop_women_public_assistance", "prop_women_labor_force", 
                      "prop_birth_last_12_months", "prop_women_HS_grad", 
                      "prop_women_college_grad", "log_occupied_housing", "log_housing_violation", 
                      "log_violent_crime", "log_nonviolent_crime")
var_individual <- c("Hispanic","Black", "Asian", "multiple_birth", "age")


p1 <- length(var_individual)
p2 <- length(var_neighborhood)
col_to_be_scaled <- p1:(p1+p2)

year <- 8; add_str <- ""; add_str2 <- "_YEAR8"
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
  assign(ifelse(output_string == "PRETERM","data_CAR_PRE","data_CAR_STILL"),data_CAR)
}

var1 <- data_CAR_PRE$alpha
var2 <- data_CAR_STILL$alpha

# Run the following twice, once for each outcome, by chaning PRE_bool to TRUE and FALSE.
PRE_bool <- TRUE # FALSE # TRUE
shift_bool <- FALSE # shift (i.e. calibrating probabilities), is needed only for SMOTE or reweight

if(PRE_bool){
  var <- var1; group_ind <- match(data_CAR_PRE$tractID, tracts@data$GEOID10)
} else {
  var <- var2; group_ind <- match(data_CAR_STILL$tractID, tracts@data$GEOID10)
}
vars <- c(var)


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

title_string <- ifelse(PRE_bool, "Preterm birth", "Stillbirth")
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

file_name <- "new_alpha"
file_name <- paste0(file_name, "_", ifelse(PRE_bool, "PRE", "STILL"),"_")
file_name <- paste0(file_name, "noSMOTE_CAR.png")
ggsave(file_name, p, scale = 1, width = 8, height = 8, units = 'in')


