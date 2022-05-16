################################################################################
######## Simulation study for Un-replicated spatial-block GEV approach #########
######## Daniela Cuba
######## April 29th, 2022
######## Functions to run simulations on the clusters
################################################################################

##### THIS IS A SECOND EXPLORATION WHERE A TRUE RESOLUTION IS SET AND THE MAXIMA WITHIN IS 
##### REPLACED WITH GEV VALUES

#### Auxiliary functions
#### Functions
# find.offset
# This function explores possible offsets of a grid resolution. It does not 
# perform offset selection. 
# Parameters: data.sp- input data as spatial points dataframe
#             resolution - the resolution considered - this must be in relation to the extent
#             off.ds - offset delta. off.ds = delta_x, delta_y

find.offset<-function(data.sp,resolution, off.ds){
  # Define extent of the data
  d.extent<-extent(data.sp)
  # Number of cells in the x and y directions (+ 1 to ensure grid is big enough even when shifting it)
  x.cells<-ceiling((d.extent[2]-d.extent[1])/resolution)+1
  y.cells<-ceiling((d.extent[4]-d.extent[3])/resolution)+1
  # Number of offsets considered in each direction
  all.offsets<- resolution/off.ds 

  #### Storage
  extent_offsets<-list()
  # maxima_offsets<-list()
  mind_offsets<-list()
  # freqs_offsets<-list()
  gev_offsets<-list()
  tests_offsets<-list()

  #### Grid exploration
  for(i in 0:all.offsets){ # shifts in the x-direction
    for(j in 0:all.offsets){ # shifts in the y-direction
      
      # Name offset
      offset_name<-paste0("x",i,".y",j)
      
      ### Build new extent
      # New extent needs to be built providing ample room near the origin
      # to allow the shifts in the grid selection
      loop_xmin<-loop_xmax-x.cells*resolution
      loop_xmax<-d.extent[2]+i*off.ds
      loop_ymax<-d.extent[4]+i*off.ds
      loop_ymin<-loop_ymax-y.cells*resolution
      loop_extent<-c(loop_xmin,loop_xmax,loop_ymin,loop_ymax)
      #extent_offsets[[offset_name]]<-loop_extent
    
      ### Produce the grid at that offset
      # Make raster
      r <- raster(ext = extent(loop_extent), res=c(resolution,resolution))
      values(r)<-1:length(r) # number grid cells
      
      # Extract Maxima
      data.sp$rcell<-raster::extract(r,data.sp)
      # freqs_offsets[[offset_name]]<-table(data.sp$rcell)
      data.max<-as.data.table(data.sp)
      data.max<-data.max[data.max[,.I[which.max(Sim)],by=rcell]$V1]
      
      ### Fit GEV
      gev_offsets[[offset_name]]<-gev.fit(data.max$Sim,show=F)
      # browser()
      # ks.test_offsets[[offset_name]]<-ks.test(as.vector(data.max$Sim),"pgev",
      #                                         gev_offsets[[offset_name]]$mle[1],
      #                                         gev_offsets[[offset_name]]$mle[2],
      #                                         gev_offsets[[offset_name]]$mle[3])
      
      ### Perform AD test - this can break sometimes so it's best in a tryCatch
      tryCatch(tests_offsets[[offset_name]]<-gnfit2(dat=as.vector(data.max$Sim),
                                                    dist="gev",
                                                    pr=gev_offsets[[offset_name]]$mle),
               error=function(cond) {return(NA)})
      
      ### Find minimum distance between maxima
      # max.sp<-SpatialPointsDataFrame(coords=as.data.frame(data.max[,3:4]),
      #                                data=as.data.frame(data.max[,1:2]))
      # maxima_offsets[[offset_name]]<-max.sp
      # max_mat<-as.matrix(dist(as.data.frame(max.sp@coords)))
      
      # Save minimum distance
      # mind_offsets[[offset_name]]<-min(as.vector(apply(max_mat,2,function(x) min(x[which(x>0)]))))
      # print(paste0("Completed j=",j," out of ",all.offsets))
    }
    # print(paste0("Completed i=",i," out of ",all.offsets))
  }
  
  # Save output
  return(list(#extents=extent_offsets,
    # maxima=maxima_offsets,
    mins=mind_offsets,
    # freqs=freqs_offsets,
    gev=gev_offsets,
    gof.tests=tests_offsets))
}

# Alteration of gnfit function to fit AD test from gnFit package
# Does away with plot and messages
# gnfit2 - fit Alderson-Darling test without the annoying plots and messages

gnfit2<-function (dat, dist, df = NULL, pr = NULL, threshold = NULL) 
{
  dat <- as.numeric(dat)
  x <- NULL
  z <- list()
  op <- par(mfrow = c(1, 2))
  if (is.null(pr)) {
    if (dist == "gev" | dist == "gpd") 
      stop("Enter Parameters!")
    else if (dist == "t") {
      if (df > 2) {
        loc <- mean(dat)
        sc <- sqrt((df - 2) * var(dat)/df)
        xdat <- (dat - loc)/sc
        prob <- pt(xdat, df)
        # qqplot(qt(ppoints(500), df), xdat, main = paste0("Q-Q plot for ", 
        #                                                  dist, ". distribution-", "DF:", 
        #                                                  df), xlab = "", ylab = "")
        # qqline(xdat, distribution = function(p) qt(p, 
        #                                            df), probs = c(0.1, 0.6), col = "blue")
        # hist(xdat, breaks = 20, prob = TRUE, xlab = "x-normalized value", 
        #      main = paste0(dist, ". distribution curve over histogram"))
        # curve(dt(x, df), col = "blue", add = TRUE, 
        #       yaxt = "n")
        # points(xdat, rep(0, length(xdat)))
      }
      else stop("DF must be > 2")
    }
    else if (dist == "gum") {
      sc <- sqrt(6 * var(dat)/pi^2)
      loc <- mean(dat) - 0.577 * sc
      pr <- c(loc, sc, 0)
      prob <- gevf(pr, dat)
      # gev.qq(pr, dat)
      # gev.his(pr, dat)
    }
    else {
      loc <- mean(dat)
      ifelse(dist == "norm", sc <- sd(dat), NA)
      ifelse(dist == "laplace", sc <- sqrt(var(dat)/2), 
             NA)
      ifelse(dist == "logis", sc <- sqrt(3 * var(dat)/pi^2), 
             NA)
      prob <- get(paste0("p", dist))(dat, loc, sc)
      # qqplot(get(paste0("q", dist))(ppoints(500), 
      #                               loc, sc), dat, main = paste0("Q-Q plot for ", 
      #                                                            dist, ". distribution"), xlab = "", 
      #        ylab = "")
      # qqline(dat, distribution = function(p) get(paste0("q", 
      #                                                   dist))(p, loc, sc), probs = c(0.1, 0.6), col = "blue")
      # hist(dat, breaks = 15, prob = TRUE, xlab = "x-normalized value", 
      #      main = paste0(dist, ". distribution curve over histogram"))
      # curve(get(paste0("d", dist))(x, loc, sc), col = "blue", 
      #       add = TRUE, yaxt = "n")
      # points(dat, rep(0, length(dat)))
    }
  }
  else {
    if (dist == "gev") {
      prob <- gevf(pr, dat)
      # gev.qq(pr, dat)
      # gev.his(pr, dat)
    }
    else if (dist == "gum") {
      pr[3] <- 0
      prob <- gevf(pr, dat)
      # gev.qq(pr, dat)
      # gev.his(pr, dat)
    }
    else if (dist == "gpd") 
      if (!is.null(threshold)) {
        u <- threshold
        dat <- dat[dat > u]
        prob <- gpdf(pr, u, dat)
        # gpd.qq(pr, u, dat)
        # gpd.his(pr, u, dat)
      }
    else stop("threshold is missing!")
  }
  n <- length(dat)
  k <- seq(1:n)
  qnor <- qnorm(sort(prob))
  pnor <- pnorm((qnor - mean(qnor))/sd(qnor))
  w <- round((sum((pnor - (2 * k - 1)/(2 * n))^2) + 1/(12 * 
                                                         n)) * (1 + 0.5/n), 4)
  if (w < 0.0275) {
    pval <- 1 - exp(-13.953 + 775.5 * w - 12542.61 * w^2)
  }
  else if (w < 0.051) {
    pval <- 1 - exp(-5.903 + 179.546 * w - 1515.29 * w^2)
  }
  else if (w < 0.092) {
    pval <- exp(0.886 - 31.62 * w + 10.897 * w^2)
  }
  else if (w < 1.1) {
    pval <- exp(1.111 - 34.242 * w + 12.832 * w^2)
  }
  else {
    warning("p-value is smaller than 7.37e-10")
  }
  z$Wpval <- pval
  A <- (-n - sum((2 * k - 1) * log(pnor) + (2 * n + 1 - 2 * 
                                              k) * log(1 - pnor))/n) * (1 + 0.75/n + 2.25/n^2)
  A <- round((1 + 0.75/n + 2.25/n^2) * A, 4)
  if (A < 0.2) {
    pval <- 1 - exp(-13.436 + 101.14 * A - 223.73 * A^2)
  }
  else if (A < 0.34) {
    pval <- 1 - exp(-8.318 + 42.796 * A - 59.938 * A^2)
  }
  else if (A < 0.6) {
    pval <- exp(0.9177 - 4.279 * A - 1.38 * A^2)
  }
  else if (A < 10) {
    pval <- exp(1.2937 - 5.709 * A + 0.0186 * A^2)
  }
  else {
    warning("p-value is smaller than 7.37e-10")
  }
  z$Apval <- pval
  z$Cram <- w
  z$Ander <- A
  # message("Test of Hypothesis for ", dist, " distribution")
  # message("Cramer-von Misses Statistics:  ", z$Cram, 
  #         "   P-Value:  ", round(z$Wpval, 5))
  # message("Anderson-Darling Statistics:   ", z$Ander, 
  #         "   P-Value:  ", round(z$Apval, 5))
  class(z) <- "gnfit"
  invisible(z)
}


# Find intersect function to identify the beginning of parameter stability
# Paramters: x is a 2-column df or matrix of 95% ci
find.intersect<-function(x){
  nrow.x<-nrow(x)
  keep.trac<-data.frame(matrix(NA,ncol=4,nrow=24))
  names(keep.trac)<-c("Change.min","Change.max","Min","Max")
  for(m in 1:(nrow(x)-1)){
    if(m==1){
      current.max<-x[nrow.x,2] #
      current.min<-x[nrow.x,1]
    }
    if(x[nrow.x-m,2]<current.max & x[nrow.x-m,2]>current.min){
      current.max<-x[nrow.x-m,2]
      hey.max<-T
    }else{
      hey.max<-F
    }
    if(x[nrow.x-m,1]>current.min & x[nrow.x-m,1]<current.max){
      current.min<-x[nrow.x-m,1]
      hey.min<-T
    }else{
      hey.min<-F
    }
    
    if((x[nrow.x-m,1]>current.max) | (x[nrow.x-m,2]<current.min)){
      is.in<-F
    }else{
      is.in<-T
    }
    keep.trac$Change.min[m]<-hey.min
    keep.trac$Change.max[m]<-hey.max
    keep.trac$Min[m]<-current.min
    keep.trac$Max[m]<-current.max
    keep.trac$Overlap[m]<-is.in
    # maxes[m]<-current.max
    # mines[m]<-current.min
  }
  keep.trac$Index<-(nrow.x-1):1
  return(keep.trac)
} 

# Parameter stability conversions functions
mu.star<-function(mu_j,j,xi,sigma_j){
  mu_j-sigma_j*(j^(-xi))*(((j^xi)-1)/xi)
}
sigma.star<-function(sigma_j,j,xi){
  sigma_j*(j^(-xi))
}

var.mustar<-function(j,xi,sigma_j,cov.mu){
  # browser()
  nabla.u<-c(1,(-j^(-xi))*((j^xi)-1)/xi,
             sigma_j*((1-(j^(-xi))*(xi*log(j)+1))/(xi^2)))
  return(t(nabla.u )%*% cov.mu %*% nabla.u)
}
var.sigmastar<-function(j,sigma_j,xi,cov.sigma){
  nabla.sigma<-c(0,j^(-xi),
                 (-sigma_j)*(j^(-xi))*log(j))
  return(t(nabla.sigma)%*% cov.sigma %*% nabla.sigma)
}


# Simulate gp-gev data
# sim.gp.gev2 - simulation gaussian process with gev replacements
# Parameters: n - total number of observations
#             ngev - number of locations replaced with gev
#             mat.nu - smoothing parameter of matern function
#             mat.range - range of matern function
#             gev.location - location of gev distribution
#             gev.scale - scale of gev distribution
#             gev.shape - shape of gev distribution
#             boxbox - bounding box of region of interest (four row matrix, with x coordinates in column 1 and y coords in column 2)
#             nsim - number of simulations
## Changes implement - the function has a "true" resolution value and
## changes the block maxima for gev values
## parameter: true.res - true resolution - either one number or a vector of lenth nsim
simulation.gp.gev<-function(n=1000,ngev,mat.nu=1,mat.range,
                             gev.location=6,gev.scale,gev.shape,
                             boxbox,nsim=10,true.res){
  # storage
  block.sel.loop<-list()
  # sim.intersect.table<-data.frame(matrix(NA,ncol=3,nrow=nsim))
  # par.intersect.table<-data.frame(matrix(NA,ncol=3,nrow=nsim))
  # Prep true.res parameter
  if(length(true.res)==1){
    true.res<-rep(true.res,
                  length.out=nsim)
  } else if(length(true.res)>1){
    true.res<-true.res
  } 
  
  # Loop to produce realizations
  for(m in 1:nsim){
    ### Produce Simulated GP
    # Simulated gp onto coordinates
    # print(paste0("Simulation n=",m))
    # browser()
    sim.gp.list<-geoR::grf(n=n,grid="irreg",xlims=range(boxbox[,1]),ylims=range(boxbox[,2]),
                     cov.model="matern", cov.pars=c(1,mat.range), kappa=mat.nu, nugget=0, 
                     lambda = 1,mean=rep(0,n), RF=F)
    print(paste0("Simulation Number ",m," of ",nsim, " at ",Sys.time()))
    sim.gp.data<-sim.gp.list$data
    sim.gp.coords<-sim.gp.list$coords
    
    
    ###### Replace block-maxima from true resolution using a GEV
    # Divide data into blocks
    # browser()
    sim.gp.gev.sp<-SpatialPointsDataFrame(coords=sim.gp.coords,
                                          data=data.frame(Sim=sim.gp.data))
    new.max<-true.res[m]*((100%/%true.res[m])+1) # +1 is just to ensure the grid is large enough to cover the data
    new.ext<-c(0,new.max,
               0,new.max)
    r.sim<-raster(ext=extent(new.ext),
                  res=true.res[m])
    values(r.sim)<-1:length(r.sim) # number grid
    sim.gp.gev.sp$rcell<-raster::extract(r.sim,sim.gp.gev.sp) # divide points
    
    # Simulate from gev distribution
    sim.gev.data<-extRemes::revd(length(r.sim),
                                 loc=gev.location,
                                 scale=gev.scale,
                                 shape=gev.shape) 
    
    # Replace maxima
    max.sim<-as.data.table(sim.gp.gev.sp@data)
    max.sim[max.sim[,.I[which.max(Sim)],by=rcell]$V1]$Sim<-sim.gev.data[sample(1:length(r.sim),
                                                                               size=length(max.sim[,.I[which.max(Sim)],by=rcell]$V1),
                                                                               replace=F)]

    #### Perform block selection
    sim.gp.gev.sp<-SpatialPointsDataFrame(coords=sim.gp.coords,
                                          data=data.frame(Sim=max.sim$Sim))
    
    # Candidate Resolutions
    res.add<-(max(boxbox[,1])-min(boxbox[,1]))/100 # increment of change
    resolutions<-seq(from=min(boxbox[,1])+res.add,
                     to=max(boxbox[,1])/4,
                     by=res.add)
    
    # Fit perform block selection ---
    #### This bit is particularly inefficient as it stores loads of unnecessary data. 
    #### It was originally coded like this to get more detailed information on the 
    #### offset exploration. Easy to make better. 
    # Go through all resolutions
    block.selection<-list()
    for(z in 1:length(resolutions)){
      loop.name<-paste0("Resolution",resolutions[z])
      block.selection[[loop.name]]<-find.offset(data.sp=sim.gp.gev.sp,
                                                resolution=resolutions[z], off.ds=0.25)
    }
    # browser()
    
    ####################################################################
    #########   Perform offset selection for each resolution   #########
    stable.mu<-list()
    stable.sigma<-list()
    sim.shape.list<-list()
    
    # Obtain mles
    sim.locs.all<-lapply(block.selection, function(x){ 
      as.vector(unlist(lapply(x$gev,function(y){ y$mle[1]})))})
    sim.scale.all<-lapply(block.selection, function(x){ 
      as.vector(unlist(lapply(x$gev,function(y){ y$mle[2]})))})
    sim.shape.all<-lapply(block.selection, function(x){
      as.vector(unlist(lapply(x$gev,function(y){ y$mle[3]})))})

    # Old code to find parameter stability using all offset information
    # ### Location
    # # Convert to mu.star
    # stable.mu<-lapply(1:length(sim.locs.all),function(x){
    #   data.frame(mu.star=mu.star(mu_j=sim.locs.all[[x]],
    #                              j=rep(x,length(sim.locs.all[[x]])),
    #                              xi=sim.shape.all[[x]],
    #                              sigma_j=sim.scale.all[[x]]),
    #              resolution=rep(x,length(sim.locs.all[[x]])))
    # })%>% rbindlist %>% as.data.frame
    # # Obtain interquartile range of all offsets per resolution
    # ranges.mu.min<-aggregate(mu.star~resolution,function(x) {quantile(x,.25)},data=stable.mu)
    # ranges.mu.max<-aggregate(mu.star~resolution,function(x) {quantile(x,.75)},data=stable.mu)
    # mu.intersects<-find.intersect(data.frame(ranges.mu.min[,2], ranges.mu.max[,2]))
    # 
    # #### Scale
    # stable.sigma<-lapply(1:length(sim.scale.all),function(x){
    #   data.frame(sigma.star=sigma.star(sigma_j=sim.scale.all[[x]],
    #                                    j=rep(x,length(sim.locs.all[[x]])),
    #                                    xi=sim.shape.all[[x]]),
    #              resolution=rep(x,length(sim.locs.all[[x]])))
    # })%>% rbindlist %>% as.data.frame
    # ranges.sigma.min<-aggregate(sigma.star~resolution,function(x) {quantile(x,.25)},data=stable.sigma)
    # ranges.sigma.max<-aggregate(sigma.star~resolution,function(x) {quantile(x,.75)},data=stable.sigma)
    # sigma.intersects<-find.intersect(data.frame(ranges.sigma.min[,2],
    #                                             ranges.sigma.max[,2]))
    # 
    # 
    # ### Shape
    # stable.shape<-lapply(1:length(sim.shape.all),function(x){
    #   data.frame(shape=sim.shape.all[[x]],
    #              resolution=rep(x,length(sim.shape.all[[x]])))
    # })%>% rbindlist %>% as.data.frame
    # 
    # ranges.shape.min<-aggregate(shape~resolution,function(x) {quantile(x,.25)},data=stable.shape)
    # ranges.shape.max<-aggregate(shape~resolution,function(x) {quantile(x,.75)},data=stable.shape)
    # shape.intersects<-find.intersect(data.frame(ranges.shape.min[,2],
    #                                             ranges.shape.max[,2]))
    # 
    # sim.intersect.table<-c(mu.intersects[max(which(mu.intersects$Overlap)),6],
    #                        sigma.intersects[max(which(sigma.intersects$Overlap)),6],
    #                        shape.intersects[max(which(shape.intersects$Overlap)),6])
    # all.intersects3[[i]]<-list(mu=mu.intersects,
    #                            sigma=sigma.intersects,
    #                            shape=shape.intersects)
    
    ######################## Best gev for each resolution only ######################
    
    # Select best offset by ADTest statistic
    best.gevs<-lapply(block.selection,
                      function(x){
                        best.offset<-lapply(x$gev,
                                            function(y) y$gof.tests$Ander) %>% unlist %>% as.vector %>% which.min
                        return(x$gev[[best.offset]])})
    
    # Extract and convert all parameters to stable parameters
    stable.pars<-lapply(1:length(sim.shape.all),function(x){
      best.offset<-lapply(block.selection[[x]]$gev,
                          function(y) y$nllh) %>% unlist %>% as.vector %>% which.min
      stable.par_v<-c(mu.star(mu_j=sim.locs.all[[x]],
                              j=rep(x,length(sim.locs.all[[x]])),
                              xi=sim.shape.all[[x]],
                              sigma_j=sim.scale.all[[x]])[best.offset],
                      sigma.star(sigma_j=sim.scale.all[[x]],
                                 j=rep(x,length(sim.locs.all[[x]])),
                                 xi=sim.shape.all[[x]])[best.offset],
                      sim.shape.all[[x]][best.offset])
      stable.par_se<-c(var.mustar(j=x,xi=block.selection[[x]]$gev[[best.offset]]$mle[3],
                                  cov.mu=block.selection[[x]]$gev[[best.offset]]$cov,
                                  sigma_j=block.selection[[x]]$gev[[best.offset]]$mle[2]) %>% sqrt,
                       var.sigmastar(j=x,xi=block.selection[[x]]$gev[[best.offset]]$mle[3],
                                     cov.sigma=block.selection[[x]]$gev[[best.offset]]$cov,
                                     sigma_j=block.selection[[x]]$gev[[best.offset]]$mle[2]) %>% sqrt,
                       block.selection[[x]]$gev[[best.offset]]$se[3])
      return(list(stable.parameters=stable.par_v,
                  stable.se=stable.par_se))
    })
    
    # Identify parameter stability resolution using modified parameter CIs
    est.mu<-lapply(stable.pars,function(x) x$stable.parameters[1])%>%unlist%>%as.vector
    se.mu<-lapply(stable.pars,function(x) x$stable.se[1]) %>% unlist %>%as.vector
    
    intersect.best.stable.pars.mu<-find.intersect(data.frame(est.mu-1.96*se.mu,
                                                             est.mu+1.96*se.mu))
    est.sigma<-lapply(stable.pars,function(x) x$stable.parameters[2])%>%unlist%>%as.vector
    se.sigma<-lapply(stable.pars,function(x) x$stable.se[2]) %>% unlist %>%as.vector
    intersect.best.stable.pars.sigma<-find.intersect(data.frame(est.sigma-1.96*se.sigma,
                                                                est.sigma+1.96*se.sigma))
    
    est.shape<-lapply(stable.pars,function(x) x$stable.parameters[3])%>%unlist%>%as.vector
    se.shape<-lapply(stable.pars,function(x) x$stable.se[3]) %>% unlist %>%as.vector
    intersect.best.stable.pars.shape<-find.intersect(data.frame(est.shape-1.96*se.shape,
                                                                est.shape+1.96*se.shape))
    par.intersect.table<-c(intersect.best.stable.pars.mu[max(which(intersect.best.stable.pars.mu$Overlap)),6],
                           intersect.best.stable.pars.sigma[max(which(intersect.best.stable.pars.sigma$Overlap)),6],
                           intersect.best.stable.pars.shape[max(which(intersect.best.stable.pars.shape$Overlap)),6])
    
    ####################################################################
    
    
    # Save output for each simulation
    block.sel.loop[[paste0("sim",m)]]<-list(Sim.GP=sim.gp.list,
                                            Sim.Gev=sim.gev.data,
                                            Sim.Mix=sim.gp.gev.sp@data$Sim,
                                            Stable.pars=stable.pars,
                                            Stable.par.intersects.all=sim.intersect.table,
                                            Stable.par.intersects.best=par.intersect.table,
                                            Block.res.sel=best.gevs,
                                            parameters=list(gev.scale=gev.scale,
                                                            gev.shape=gev.shape,
                                                            gev.loc=gev.location,
                                                            mat.range=mat.range))
    
  }
  
  # Return correct output
  # save(block.sel.loop,file=save.as)
  return(block.sel.loop)
}

