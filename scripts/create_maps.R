# setwd("../cecis_scripts/")
library(sp)
library(spdep)
library(rgeos)
library(RColorBrewer)
library(ggplot2)
library(ggmap)
library(dplyr)
library(maptools)

# load("../../../data/shapefile/phillytracts")
tmp = load("../../../data/phillygooglemap.rdata")
wdstr <- "results/results_feb8_9/"
burnin <- 500
thin <- 5
mcmc_niter <- 2001
index <- seq(burnin, mcmc_niter, by = thin)

SMOTE_bool <- FALSE

if(SMOTE_bool){
  tmp <- load(paste0(wdstr,"output_alldata_nogamma_rho_SMOTE_nointeractions_newcov_PRETERM_CAR_correct.RData"))
} else {
  tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_newcov_PRETERM_CAR_correct.RData"))
}

data_CAR_PRE <- data.frame(tractID = tractID_unique_list[[1]],  
                           phat_all = mean_phat_all_list[[1]], phat2_all = mean_phat2_all_list[[1]],
                           phat2shift_all = mean_phat2shift_all_list[[1]], phat_neig = phat_neig_list[[1]],
                           phatshift_neigh = phat_neigshift_list[[1]],
                           alpha = rowMeans(output_list[[1]]$alpha[,index]))
rm(list = tmp)

if(SMOTE_bool){
  tmp <- load(paste0(wdstr,"output_alldata_nogamma_rho_SMOTE_nointeractions_newcov_STILLBIRTH_CAR_correct.RData"))
} else {
  tmp <- load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_newcov_STILLBIRTH_CAR_correct.RData"))
}
data_CAR_STILL <- data.frame(tractID = tractID_unique_list[[1]],  
                             phat_all = mean_phat_all_list[[1]], phat2_all = mean_phat2_all_list[[1]],
                             phat2shift_all = mean_phat2shift_all_list[[1]], phat_neig = phat_neig_list[[1]],
                             phatshift_neigh = phat_neigshift_list[[1]],
                             alpha = rowMeans(output_list[[1]]$alpha[,index]))
rm(list = tmp)



# if(PRE_bool){
#   tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_newcov_PRETERM_CAR_correct.RData"))
# } else {
#   tmp <-load(paste0(wdstr,"output_alldata_nogamma_rho_nointeractions_newcov_STILLBIRTH_CAR_correct.RData"))
# }
# 
# data_CAR_PRE_noSMOTE <- data.frame(tractID = tractID_unique_list[[1]],  
#                                    phat_all = mean_phat_all_list[[1]], phat2_all = mean_phat2_all_list[[1]],
#                                    phat2shift_all = mean_phat2shift_all_list[[1]], phat_neig = phat_neig_list[[1]],
#                                    phatshift_neigh = phat_neigshift_list[[1]])
# rm(list = tmp)
# 
# if(PRE_bool){
#   tmp <- load(paste0(wdstr,"output_alldata_nogamma_rho_SMOTE_nointeractions_newcov_PRETERM_CAR_correct.RData"))
# } else {
#   tmp <- load(paste0(wdstr,"output_alldata_nogamma_rho_SMOTE_nointeractions_newcov_STILLBIRTH_CAR_correct.RData"))
# }
# data_CAR_PRE_SMOTE <- data.frame(tractID = tractID_unique_list[[1]],  
#                                    phat_all = mean_phat_all_list[[1]], phat2_all = mean_phat2_all_list[[1]],
#                                    phat2shift_all = mean_phat2shift_all_list[[1]], phat_neig = phat_neig_list[[1]],
#                                    phatshift_neigh = phat_neigshift_list[[1]])
# rm(list = tmp)


PRE_bool <- TRUE

alpha_bool <- TRUE
phat_bool <- FALSE
shift_bool <- FALSE

if(alpha_bool){
  var1 <- data_CAR_PRE$alpha
  var2 <- data_CAR_STILL$alpha
} else {
  if(phat_bool){
    var1 <- data_CAR_PRE$phat_all
    var2 <- data_CAR_STILL$phat_all
  } else {
    if(shift_bool){
      var1 <- data_CAR_PRE$phat2shift_all
      var2 <- data_CAR_STILL$phat2shift_all
    } else {
      var1 <- data_CAR_PRE$phat2_all
      var2 <- data_CAR_STILL$phat2_all
    }
  }
}

if(PRE_bool){
  var <- var1; group_ind <- match(data_CAR_PRE$tractID, tracts@data$GEOID10)
} else {
  var <- var2; group_ind <- match(data_CAR_STILL$tractID, tracts@data$GEOID10)
}
# vars <- c(var1, var2)
vars <- c(var) # for now we ignore the comparison with SMOTE

# if(phat_bool){
#   var1 <- data_CAR_PRE_SMOTE$phat_all
#   var2 <- data_CAR_PRE_noSMOTE$phat_all
# } else {
#   if(shift_bool){
#     var1 <- data_CAR_PRE_SMOTE$phat2shift_all
#     var2 <- data_CAR_PRE_noSMOTE$phat2shift_all
#   } else {
#     var1 <- data_CAR_PRE_SMOTE$phat2_all
#     var2 <- data_CAR_PRE_noSMOTE$phat2_all
#   }
# }
# if(SMOTE_bool){
#   var <- var1; group_ind <- match(data_CAR_PRE_SMOTE$tractID, tracts@data$GEOID10)
# } else {
#   var <- var2; group_ind <- match(data_CAR_PRE_noSMOTE$tractID, tracts@data$GEOID10)
# }
# # vars <- c(var1, var2)
# vars <- c(var) # for now we ignore the comparison with SMOTE




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
if(alpha_bool){
  p <- p + scale_fill_distiller(palette = palette, limits = c(MI,MA), values = c(0,z_adj,1), name =expression(alpha[i]))
} else {
  p <- p + scale_fill_distiller(palette = palette, limits = c(MI,MA), values = c(0,z_adj,1), name =expression(hat(p)[i]))
}
p <- p + theme(plot.title = element_text(hjust = 0.5, vjust = -1,size=16),
        panel.border =  element_rect(colour = 'transparent', fill = 'transparent'),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, 'mm'), axis.text = element_blank(),
        panel.grid = element_blank(), axis.title = element_blank(),
        legend.key.height = unit(1,"cm"),
        legend.text=element_text(size=12),
        legend.title = element_text(size=15),
        legend.justification = c(1,0.035), legend.position = c(1,0.035),
        legend.background = element_rect(fill= alpha("white",0.6)))
p <- p + scale_bar(lon = -75.17, lat = 39.825, 
                         distance_lon = 5, distance_lat = 0.5, distance_legend = 1.1, 
                         dist_unit = "km", orientation = FALSE)
# p
file_name <- ifelse(alpha_bool, "alpha", ifelse(phat_bool, "phat", ifelse(shift_bool, "phat2shift", "phat2")))
file_name <- paste0(file_name, "_", ifelse(PRE_bool, "PRE", "STILL"),"_")
file_name <- paste0(file_name, ifelse(SMOTE_bool, "SMOTE", "noSMOTE"),"_")
file_name <- paste0(file_name, "CAR.png")
ggsave(file_name, p, scale = 1, width = 8, height = 8, units = 'in')

