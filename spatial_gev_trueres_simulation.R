#########################################################################
################### Univariate Spatial GEV - simulations ################
### Author: Daniela Cuba
### Date: April 29th, 2022
### Script to run in the cluster for simulations
#########################################################################


# Load libraries
libs<-c("sp","evd","fields","data.table","tidyverse","rgdal","ismev",
        "raster","gnFit")
lapply(libs,require,character.only=T)

# Load functions
source("./code/spatial_gev_trueres_cluster_functions.R")

# Run code
library(parallel)
B<-1000
boxbox<-matrix(c(0,0,0,100,100,0,100,100),byrow=T,ncol=2)
save.progress<-list()
for(B in 1:(1000%/%10)){
  sim.spgev<-mclapply(1:10, function(x){
    sim.spgev_y<-simulation.gp.gev(n=1000, 
                                   ngev=250,
                                   mat.range=10,
                                   gev.location=6,
                                   gev.scale=1,
                                   gev.shape=0.5,
                                   boxbox=boxbox,
                                   nsim=1,
                                   true.res=8)
    },mc.cores=10)
  save.progress<-c(save.progress,sim.spgev)
  save(save.progress,file=paste0("./output/sim_output_parset1.RData"))
  print(paste0(B*10," of 1000 at ", Sys.time()))
}
# save(save.progress,file="./output/sim_output_parset1.RData")
