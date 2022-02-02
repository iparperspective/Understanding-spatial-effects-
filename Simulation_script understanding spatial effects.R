



#######################################
### simulate through point-process ####
#######################################

require(knitr)
require(INLA)
require(rgdal)
require(RandomFields)
require(lattice)
require(scales)
require(gridExtra)
library(grid)
require(ggplot2)
require(overlapping)
require(scoringRules)
require(raster)
require(tidyr)
require(xtable)
require(mgcv)

RFgetModelNames(type="positive definite", domain="single variable",
                iso="isotropic") 

#INLA:::inla.dynload.workaround()
inla.setOption(scale.model.default = TRUE)

### SSIM funcitons
{
  # Appendix C #
  
  
  library(maptools)
  library(geepack)
  library(splines)
  library(SpatialTools)
  
  
  #*******************************************************************************
  #SSIM - R implementation of the Structural Similarity Index using raster package
  #*******************************************************************************
  
  library(raster)  ### Requires raster >= 2.3-12
  
  #Function for iterative edge correction via map reflection (and averaging)
  edge.cor.ref <- function(ras,w){
    iter.edge.cor.ref <- function(ind, ras, cel){
      #reflect along vertical/horizontal edges
      i <- cel*2-1
      col. <- colFromCell(ras,ind)
      row. <- rowFromCell(ras,ind)
      cols <- c(col.-i, col.,   col.+i, col.)
      rows <- c(row.,   row.-i, row.,   row.+i)
      celly <- cellFromRowCol(ras,rows,cols)
      subr <- ras[celly]
      est <- mean(subr,na.rm=TRUE)
      if (!is.nan(est)){
        return(est)
      } else {
        #Reflect along diagonal if it is only corner connected.
        cols <- c(col.-i, col.-i, col.+i, col.+i)
        rows <- c(row.-i, row.+i, row.+i, row.-i)
        celly <- cellFromRowCol(ras,rows,cols)
        subr <- ras[celly]
        est <- mean(subr,na.rm=TRUE)
        if (!is.nan(est)){
          return(est)
        } else {return(NA)}
      }
    }
    #--------------
    #Assume NA values (edge padding or donut holes) exist for edge correction.
    for (cel in 1:w){
      temp <- boundaries(ras,type='outer')
      loc <- Which(temp==1,cells=T)
      for (ind in loc){
        ras[ind] <- iter.edge.cor.ref(ind,ras,cel)
      }
    }
    return(ras)
  }
  
  
  
  #Gaussian filter weights matrix
  filter.g <- function(w,sigma){
    f.g <- function(x,y,sigma){ (1/(2*pi*sigma^2))*exp(-(x^2+y^2)/(2*sigma^2))}
    w.i <- seq(-w,w,1)
    xy <- expand.grid(x=w.i,y=w.i)
    xy$w <- f.g(xy$x,xy$y,sigma)
    w.m <- matrix(xy$w,nrow=length(w.i),byrow=T)/sum(xy$w)
    return(w.m)
  }
  
  #==========================================
  ssimMap <- function(img1, img2, w=3, sigma=1.5, gFil=FALSE, outer.edge.pad=FALSE, edge.cor=FALSE) {
    
    #Check to see if extents are equal
    img1.extent <- extent(img1)
    img2.extent <- extent(img2)
    img1.na <- Which(is.na(img1),cells=TRUE)
    if (img1.extent != img2.extent){stop('Warning: SSIM calculation aborted. The raster extents do not match.')}
    
    #set constants
    l <- max(cellStats(img1, max), cellStats(img2, max))
    globalMin <- abs(min(cellStats(img1, min), cellStats(img2, min)))
    l <- l - globalMin
    k <- c(0.01, 0.03)
    C1 <-(k[1]*l)^2
    C2 <-(k[2]*l)^2
    C3 <-C2/2
    
    #Create Null filter
    filterx <- matrix(1,ncol=w*2+1,nrow=w*2+1)/(w*2+1)^2
    if(gFil) {
      #create Gaussian filter
      filterx <- filter.g(w,sigma)
    } 
    
    #Optionally pad edges with NA's for edge correction
    if (outer.edge.pad){
      img1 <- extend(img1,2*w)
      img2 <- extend(img2,2*w)
    }
    
    #Compute iterative edge correction 'reflect'
    if (edge.cor=='reflect'){
      img1 <- edge.cor.ref(img1,w)
      img2 <- edge.cor.ref(img2,w)
    } 
    
    #get mu 
    mu1 <- raster::focal(img1, filterx)
    mu2 <- raster::focal(img2, filterx)
    
    sig1 <- abs(raster::focal(img1*img1,filterx) - mu1*mu1)^0.5
    sig2 <- abs(raster::focal(img2*img2,filterx) - mu2*mu2)^0.5
    #sig12 relates to correlation
    sig12 <- raster::focal(img1*img2, filterx) - mu1*mu2
    
    #compute components
    L <- ((2*mu1*mu2)+C1) / (mu1^2 + mu2^2 + C1)
    C <- ((2*sig1*sig2)+C2) / (sig1^2 + sig2^2 + C2)
    S <- (sig12 + C3) / (sig1 * sig2 + C3)
    #compute SSIM
    SSIM2 <- L * C * S
    
    #Compute RasterBrick
    ssim.brick <- brick(SSIM2, L, C, S)  
    ssim.brick <- raster::crop(ssim.brick,img1.extent)
    ssim.brick[img1.na] <- NA
    
    ssim.brick@data@names <- c('SSIM', 'SIM', 'SIV', 'SIP')
    return(ssim.brick)
  }
  
  
  
  
}

## spatial autocor due to proximity
# also avoid 0 probability locations
ag1=.15
ag2=.15
min.prob = .1

## define the locations:
from <- 0
to <- 10
#to <- 5
x.seq <- seq(from, to, length=31) 
y.seq <- seq(from, to, length=31)

## model covar space
model_small_cov <-   RMmatern(nu=2,var=2,scale=2)
model_med_cov <-   RMmatern(nu=3,var=3,scale=3)
model_large_cov <-   RMmatern(nu=10,var=10,scale=10)

n_sim=50
delta_mae_gam = delta_aic = delta_mae = delta_waic = within_rmse_gam = within_rmse_inla = outof_rmse_gam = outof_rmse_inla = within_mae_gam = within_mae_inla = outof_mae_gam = outof_mae_inla = matrix(NA,nrow = n_sim,ncol=9)
delta_M_0_SIP = delta_M_S_SIP = delta_M_M_SIP = delta_M_L_SIP = delta_M_SM_SIP = delta_M_SL_SIP = delta_M_ML_SIP = delta_M_SML_SIP = matrix(NA,nrow = n_sim,ncol=8)
delta_M_0_SIP_gam = delta_M_S_SIP_gam = delta_M_M_SIP_gam = delta_M_L_SIP_gam = delta_M_SM_SIP_gam = delta_M_SL_SIP_gam = delta_M_ML_SIP_gam = delta_M_SML_SIP_gam = matrix(NA,nrow = n_sim,ncol=8)
for(sim in 1:1){
  
  time0 = Sys.time()
  # set.seed(555)
  seed = sample(1:9999,1)
  set.seed(seed)
  cov_small_surface <- RFsimulate(model_small_cov, x=x.seq, y=y.seq)
  spplot(cov_small_surface)
  cov_small_surface$variable1 = rescale(cov_small_surface$variable1,to=c(-2,2))
  
  cov_med_surface <- RFsimulate(model_med_cov, x=x.seq, y=y.seq)
  spplot(cov_med_surface)
  cov_med_surface$variable1 = rescale(cov_med_surface$variable1,to=c(-2,2))
  
  large_cov_surface  <- RFsimulate(model_large_cov, x=x.seq, y=y.seq) 
  spplot(large_cov_surface)
  large_cov_surface$variable1 = rescale(large_cov_surface$variable1,to=c(-2,2))
  
  
  strength_small_cov = 1
  strength_med_cov = 1
  strength_large_cov = - 1
  
  effect_small_cov = cov_small_surface$variable1*strength_small_cov
  effect_med_cov = cov_med_surface$variable1*strength_med_cov
  effect_large_cov = large_cov_surface$variable1*strength_large_cov
  
  
  coords = cbind(coordinates(cov_small_surface)[,1],coordinates(cov_small_surface)[,2])
  min.dist=max(abs(diff(x.seq)))
  neighbors=sapply(1:nrow(coords),function(x){
    which(pointDistance(coords[x,],coords,lonlat=F)<(min.dist*2))
  })
  
  n.cells=length(x.seq)*length(y.seq)
  mean_sp1_per_loc=20
  n.sim=1000
  n.per.t.sp1=round((n.cells*mean_sp1_per_loc)/n.sim)
  #n.per.t.sp2=round((n.cells*mean_sp2_per_loc)/n.sim)
  abu_sp1=rep(0,n.cells)
  for(t in 1:n.sim){
    #########
    ### sp1
    #########
    
    ### without covariates
    prob_sp1 = min.prob + effect_small_cov + effect_med_cov +  effect_large_cov + sapply(1:nrow(coords),function(x){sum(abu_sp1[neighbors[[x]]]*ag1)})
    
    if(sum(prob_sp1)==min.prob*n.cells){prob_sp1=rep(1/n.cells,n.cells)}else{
      prob_sp1=(prob_sp1-min(prob_sp1)+min.prob);prob_sp1=prob_sp1/sum(prob_sp1) 
    }
    
    idx.1=sample(1:length(prob_sp1),n.per.t.sp1,prob=prob_sp1,replace = T)
    ind=table(idx.1)
    cells=as.numeric(names(ind))
    numbers=as.vector(ind)
    
    abu_sp1[cells]=abu_sp1[cells]+numbers
    
  }
  
  df.sim=data.frame(real=abu_sp1,
                    x=coordinates(cov_small_surface)[,1],
                    y=coordinates(cov_small_surface)[,2],
                    cov_S=cov_small_surface$variable1*strength_small_cov,
                    cov_M=cov_med_surface$variable1*strength_med_cov,
                    cov_L=large_cov_surface$variable1*strength_large_cov,
                    ef_small=effect_small_cov,
                    ef_med=effect_med_cov,
                    ef_large=effect_large_cov,
                    ef_all = effect_small_cov+effect_med_cov+effect_large_cov,
                    ef_SM = effect_small_cov+effect_med_cov,
                    ef_SL = effect_small_cov+effect_large_cov,
                    ef_ML = effect_large_cov+effect_med_cov,
                    resid = scale(abu_sp1)-scale(effect_small_cov+effect_med_cov+effect_large_cov))
  
  ea_plot = df.sim[sample(1:nrow(df.sim),100),]
  
  S_ef = ggplot(df.sim,aes(x=x,y=y,fill=cov_S)) + geom_tile() + ggtitle("S") +   
    scale_fill_gradient(low="black", high="white") +
    geom_point(data=ea_plot,aes(x=x,y=y),shape = 3,color="white") +
    theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5),
          legend.title=element_blank(), 
          legend.position = "none",
          #legend.text=element_text(size=15),
          axis.title=element_blank(),
          strip.text.x=element_text(size=20),
          axis.ticks=element_blank(), 
          axis.text=element_blank(),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5, linetype = "solid"))
  
  M_ef = ggplot(df.sim,aes(x=x,y=y,fill=cov_M)) + geom_tile() + ggtitle("M")+   
    scale_fill_gradient(low="black", high="white") +
    geom_point(data=ea_plot,aes(x=x,y=y),shape = 3,color="white") +
    theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5),
          legend.title=element_blank(), 
          legend.position = "none",
          #legend.text=element_text(size=15),
          axis.title=element_blank(),
          strip.text.x=element_text(size=20),
          axis.ticks=element_blank(), 
          axis.text=element_blank(),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5, linetype = "solid"))
  
  L_ef = ggplot(df.sim,aes(x=x,y=y,fill=cov_L)) + geom_tile() + ggtitle("L")+   
    scale_fill_gradient(low="black", high="white") +
    geom_point(data=ea_plot,aes(x=x,y=y),shape = 3,color="white") +
    theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5),
          legend.title=element_blank(), 
          legend.position = "none",
          #legend.text=element_text(size=15),
          axis.title=element_blank(),
          strip.text.x=element_text(size=20),
          axis.ticks=element_blank(), 
          axis.text=element_blank(),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5, linetype = "solid"))
  
  abundance_ef = ggplot(df.sim,aes(x=x,y=y,fill=real)) + geom_tile() + ggtitle("Real pattern")+   
    scale_fill_gradient(low="black", high="white") +
    geom_point(data=ea_plot,aes(x=x,y=y),shape = 3,color="white") +
    theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5),
          legend.title=element_blank(), 
          legend.position = "none",
          #legend.text=element_text(size=15),
          axis.title=element_blank(),
          strip.text.x=element_text(size=20),
          axis.ticks=element_blank(), 
          axis.text=element_blank(),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5, linetype = "solid"))
  
  resid_ef = ggplot(df.sim,aes(x=x,y=y,fill=resid)) + geom_tile() + ggtitle("Residual aggregation")+   
    scale_fill_gradient(low="black", high="white") +
    geom_point(data=ea_plot,aes(x=x,y=y),shape = 3,color="white") +
    theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5),
          legend.title=element_blank(), 
          legend.position = "none",
          #legend.text=element_text(size=15),
          axis.title=element_blank(),
          strip.text.x=element_text(size=20),
          axis.ticks=element_blank(), 
          axis.text=element_blank(),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5, linetype = "solid"))
  
  
  SM_ef = ggplot(df.sim) + geom_point(aes(x=x,y=y,color=ef_SM),size=5)
  SL_ef = ggplot(df.sim) + geom_point(aes(x=x,y=y,color=ef_SL),size=5)
  ML_ef = ggplot(df.sim) + geom_point(aes(x=x,y=y,color=ef_ML),size=5)
  SML_ef = ggplot(df.sim) + geom_point(aes(x=x,y=y,color=ef_all),size=5)
  
  #grid.arrange(S_ef,M_ef,L_ef,abundance_ef)
  
  
  
  #####################################################################################
  ############################## set modelling structure ############################
  #####################################################################################
  
  ############
  ### dataset
  n=100
  sample.idx=sample(1:nrow(df.sim),n)
  sample_data = df.sim[sample.idx,]
  pred_data = df.sim[-sample.idx,]
  df.pred = df.sim[,c("x","y","cov_S","cov_M","cov_L")]
  df.pred_outofsample = df.sim[-sample.idx,c("x","y","cov_S","cov_M","cov_L")]
  
  #######################
  #### modelling MGCV ####
  #######################
  {
    ##############
    ### No cov
    #############
    no_cov_form <- real ~  te(x,y,bs="tp")
    no_cov_gam = gam(no_cov_form,data = sample_data,method = "REML",family = "poisson")
    no_cov_gam_df = data.frame(x=df.sim$x,y=df.sim$y,space = predict.gam(no_cov_gam,df.pred,type = "response"))
    pred_no_cov_gam_spat = ggplot(no_cov_gam_df,aes(x=x,y=y,fill=space)) + geom_tile() 
    
    
    ##############
    ### cov_small
    ##############
    
    S_cov_form <- real ~  cov_S + te(x,y,bs="tp")
    S_cov_gam = gam(S_cov_form,data = sample_data,method = "REML",family = "poisson")
    S_cov_gam_df = data.frame(x=df.sim$x,y=df.sim$y,space = predict.gam(S_cov_gam,df.pred,type = "response"))
    pred_S_cov_gam_spat = ggplot(S_cov_gam_df,aes(x=x,y=y,fill=space)) + geom_tile() 
    
    ##############
    ### cov_medium
    ##############
    
    M_cov_form <- real ~  cov_M + te(x,y,bs="tp")
    M_cov_gam = gam(M_cov_form,data = sample_data,method = "REML",family = "poisson")
    M_cov_gam_df = data.frame(x=df.sim$x,y=df.sim$y,space = predict.gam(M_cov_gam,df.pred,type = "response"))
    pred_M_cov_gam_spat = ggplot(M_cov_gam_df,aes(x=x,y=y,fill=space)) + geom_tile() 
    
    ##############
    ### cov_large
    ##############
    
    L_cov_form <- real ~  cov_L + te(x,y,bs="tp")
    L_cov_gam = gam(L_cov_form,data = sample_data,method = "REML",family = "poisson")
    L_cov_gam_df = data.frame(x=df.sim$x,y=df.sim$y,space = predict.gam(L_cov_gam,df.pred,type = "response"))
    pred_L_cov_gam_spat = ggplot(L_cov_gam_df,aes(x=x,y=y,fill=space)) + geom_tile() 
    
    ##############
    ### SM cov
    ##############
    
    SM_cov_form <- real ~  cov_M + cov_S + te(x,y,bs="tp")
    SM_cov_gam = gam(SM_cov_form,data = sample_data,method = "REML",family = "poisson")
    SM_cov_gam_df = data.frame(x=df.sim$x,y=df.sim$y,space = predict.gam(SM_cov_gam,df.pred,type = "response"))
    pred_SM_cov_gam_spat = ggplot(SM_cov_gam_df,aes(x=x,y=y,fill=space)) + geom_tile() 
    
    ##############
    ### SL cov
    ##############
    
    SL_cov_form <- real ~  cov_L + cov_S + te(x,y,bs="tp")
    SL_cov_gam = gam(SL_cov_form,data = sample_data,method = "REML",family = "poisson")
    SL_cov_gam_df = data.frame(x=df.sim$x,y=df.sim$y,space = predict.gam(SL_cov_gam,df.pred,type = "response"))
    pred_SL_cov_gam_spat = ggplot(SL_cov_gam_df,aes(x=x,y=y,fill=space)) + geom_tile() 
    
    ##############
    ### ML cov
    ##############
    
    ML_cov_form <- real ~  cov_L + cov_M + te(x,y,bs="tp")
    ML_cov_gam = gam(ML_cov_form,data = sample_data,method = "REML",family = "poisson")
    ML_cov_gam_df = data.frame(x=df.sim$x,y=df.sim$y,space = predict.gam(ML_cov_gam,df.pred,type = "response"))
    
    ##############
    ### SML cov
    ##############
    
    SML_cov_form <- real ~  cov_L + cov_M + cov_S + te(x,y,bs="tp")
    SML_cov_gam = gam(SML_cov_form,data = sample_data,method = "REML",family = "poisson")
    SML_cov_gam_df = data.frame(x=df.sim$x,y=df.sim$y,space = predict.gam(SML_cov_gam,df.pred,type = "response"))
    pred_SML_cov_gam_spat = ggplot(SML_cov_gam_df,aes(x=x,y=y,fill=space)) + geom_tile() 
    
    ##############
    ### only cov
    ##############
    
    only_cov_form <- real ~  cov_L + cov_M + cov_S 
    only_cov_gam = gam(only_cov_form,data = sample_data,method = "REML",family = "poisson")
    only_cov_gam_df = data.frame(x=df.sim$x,y=df.sim$y,space = predict.gam(only_cov_gam,df.pred,type = "response"))
    pred_only_cov_gam_spat = ggplot(only_cov_gam_df,aes(x=x,y=y,fill=space)) + geom_tile() 
    
  }
  
  #######################
  #### modelling INLA ####
  #######################
  {
    d=cbind(sample_data$x,sample_data$y)
    #### prepare spatial effect
    ### prepare mesh
    long_dist <- max(diff(range(df.sim$x)),
                     diff(range(df.sim$y)))
    
    max.edge <- long_dist*.05
    mesh <- inla.mesh.2d(
      cbind(df.sim$x, df.sim$y), 
      max.edge = max.edge*c(1,2),
      cutoff = max.edge*0.3,
      offset = max.edge*c(0.5, 3.5))
    #plot(mesh, asp=1)
    
    idx_in_mesh = which(mesh$loc[,1]<10.1&mesh$loc[,1]>-0.1&mesh$loc[,2]<10.1&mesh$loc[,2]>-0.1)
    
    ### SPDE
    spde <- inla.spde2.pcmatern(
      mesh=mesh,
      prior.range = c(5, .5), ### p(r < 2) = 0.5
      prior.sigma = c(.5, .1)) ### p(sd > 1) = 0.1
    
    A <- inla.spde.make.A(
      mesh=mesh, 
      loc=cbind(sample_data$x, sample_data$y))
    A_pred <- inla.spde.make.A(
      mesh=mesh, 
      loc=cbind(pred_data$x, pred_data$y))
    
    ### prepare data stacks 
    data.stack  <- inla.stack(
      tag='est',
      data=list(y=sample_data$real), 
      A=list(1,A),
      effects=list(list(b0=rep(1, nrow(sample_data)),
                        cov_small = sample_data$cov_S,
                        cov_med = sample_data$cov_M,
                        cov_large = sample_data$cov_L), 
                   spat=1:spde$n.spde)) 
    
    pred.stack  <- inla.stack(
      tag='pred',
      data=list(y=rep(NA,nrow(pred_data))), 
      A=list(1,A_pred),
      effects=list(list(b0=rep(1, nrow(pred_data)),
                        cov_small = pred_data$cov_S,
                        cov_med = pred_data$cov_M,
                        cov_large = pred_data$cov_L), 
                   spat=1:spde$n.spde)) 
    
    both <- inla.stack(data.stack,pred.stack)
    
    idx_pred = inla.stack.index(both,tag="pred")$data
    idx_est = inla.stack.index(both,tag="est")$data
    
    ################################################################
    ############################## fit models ############################
    ##################################################################
    
    ##############
    ### No cov
    #############
    no_cov_form <- y ~ 0 + b0 +  
      f(spat, model=spde)   
    
    #### fit model
    no_cov_model <- inla(  no_cov_form,
                           family="poisson",
                           data=inla.stack.data(both),
                           control.compute=list(waic=T,cpo=T),
                           control.predictor=list( A=inla.stack.A(both), compute=TRUE, link=1),
                           verbose=FALSE,  #control.inla = ci,
                           num.threads = 1)
    
    no_cov_model_df = data.frame(x=mesh$loc[,1],y=mesh$loc[,2],space = no_cov_model$summary.random$spat$mean)
    no_cov_model_df=no_cov_model_df[idx_in_mesh,]
    no_cov_model_spat = ggplot(no_cov_model_df) + geom_point(aes(x=x,y=y,color = space),size=5) + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    SML_ef = SML_ef + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    #grid.arrange(no_cov_model_spat,SML_ef)
    
    ##############
    ### cov_small
    ##############
    
    S_cov_form <- y ~ 0 + b0 + cov_small + 
      f(spat, model=spde)   
    
    #### fit model
    S_cov_model <- inla(  S_cov_form,
                          family="poisson",
                          data=inla.stack.data(both),
                          control.compute=list(waic=T,cpo=T),
                          control.predictor=list( A=inla.stack.A(both), compute=TRUE, link=1),
                          verbose=FALSE,  #control.inla = ci,
                          num.threads = 1)
    
    S_cov_model$summary.fixed
    
    S_cov_model_df = data.frame(x=mesh$loc[,1],y=mesh$loc[,2],space = S_cov_model$summary.random$spat$mean)
    S_cov_model_df=S_cov_model_df[idx_in_mesh,]
    S_cov_model_spat = ggplot(S_cov_model_df) + geom_point(aes(x=x,y=y,color = space),size=5)  + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    ML_ef = ML_ef + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    
    #grid.arrange(S_cov_model_spat,ML_ef)
    
    
    ##############
    ### large_cov
    ##############
    L_cov_form <- y ~ 0 + b0 + cov_large + 
      f(spat, model=spde)   
    
    #### fit model
    L_cov_model <- inla(  L_cov_form,
                          family="poisson",
                          data=inla.stack.data(both),
                          control.compute=list(waic=T,cpo=T),
                          control.predictor=list( A=inla.stack.A(both), compute=TRUE, link=1),
                          verbose=FALSE,  #control.inla = ci,
                          num.threads = 1)
    
    L_cov_model$summary.fixed
    
    L_cov_model_df = data.frame(x=mesh$loc[,1],y=mesh$loc[,2],space = L_cov_model$summary.random$spat$mean)
    L_cov_model_df=L_cov_model_df[idx_in_mesh,]
    L_cov_model_spat = ggplot(L_cov_model_df) + geom_point(aes(x=x,y=y,color = space),size=5)+ geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    SM_ef = SM_ef + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    
    #grid.arrange(L_cov_model_spat,SM_ef)
    
    ##############
    ### M_cov
    ##############
    M_cov_form <- y ~ 0 + b0 + cov_med + 
      f(spat, model=spde)   
    
    #### fit model
    M_cov_model <- inla(  M_cov_form,
                          family="poisson",
                          data=inla.stack.data(both),
                          control.compute=list(waic=T,cpo=T),
                          control.predictor=list( A=inla.stack.A(both), compute=TRUE, link=1),
                          verbose=FALSE,  #control.inla = ci,
                          num.threads = 1)
    
    M_cov_model$summary.fixed
    
    M_cov_model_df = data.frame(x=mesh$loc[,1],y=mesh$loc[,2],space = M_cov_model$summary.random$spat$mean)
    M_cov_model_df=M_cov_model_df[idx_in_mesh,]
    M_cov_model_spat = ggplot(M_cov_model_df) + geom_point(aes(x=x,y=y,color = space),size=5)+ geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    SL_ef = S_ef + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    
    #grid.arrange(M_cov_model_spat,SL_ef)
    
    
    ##############
    ### SL cov
    ##############
    SL_cov_form <- y ~ 0 + b0 + cov_small + cov_large +
      f(spat, model=spde)   
    
    #### fit model
    SL_cov_model <- inla(  SL_cov_form,
                           family="poisson",
                           data=inla.stack.data(both),
                           control.compute=list(waic=T,cpo=T),
                           control.predictor=list( A=inla.stack.A(both), compute=TRUE, link=1),
                           verbose=FALSE,  #control.inla = ci,
                           num.threads = 1)
    
    SL_cov_model$summary.fixed
    
    SL_cov_model_df = data.frame(x=mesh$loc[,1],y=mesh$loc[,2],space = SL_cov_model$summary.random$spat$mean)
    SL_cov_model_df=SL_cov_model_df[idx_in_mesh,]
    SL_cov_model_spat = ggplot(SL_cov_model_df) + geom_point(aes(x=x,y=y,color = space),size=5) + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    M_ef = M_ef + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    
    #grid.arrange(SL_cov_model_spat,M_ef)
    
    ##############
    ### SM cov
    ##############
    SM_cov_form <- y ~ 0 + b0 + cov_small + cov_med +
      f(spat, model=spde)   
    
    #### fit model
    SM_cov_model <- inla(  SM_cov_form,
                           family="poisson",
                           data=inla.stack.data(both),
                           control.compute=list(waic=T,cpo=T),
                           control.predictor=list( A=inla.stack.A(both), compute=TRUE, link=1),
                           verbose=FALSE,  #control.inla = ci,
                           num.threads = 1)
    
    SM_cov_model$summary.fixed
    
    SM_cov_model_df = data.frame(x=mesh$loc[,1],y=mesh$loc[,2],space = SM_cov_model$summary.random$spat$mean)
    SM_cov_model_df=SM_cov_model_df[idx_in_mesh,]
    SM_cov_model_spat = ggplot(SM_cov_model_df) + geom_point(aes(x=x,y=y,color = space),size=5) + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    L_ef = L_ef + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    
    #grid.arrange(SM_cov_model_spat,L_ef)
    
    
    ##############
    ### ML cov
    ##############
    ML_cov_form <- y ~ 0 + b0  + cov_med + cov_large +
      f(spat, model=spde)   
    
    #### fit model
    ML_cov_model <- inla(  ML_cov_form,
                           family="poisson",
                           data=inla.stack.data(both),
                           control.compute=list(waic=T,cpo=T),
                           control.predictor=list( A=inla.stack.A(both), compute=TRUE, link=1),
                           verbose=FALSE,  #control.inla = ci,
                           num.threads = 1)
    
    ML_cov_model$summary.fixed
    
    ML_cov_model_df = data.frame(x=mesh$loc[,1],y=mesh$loc[,2],space = ML_cov_model$summary.random$spat$mean)
    ML_cov_model_df=ML_cov_model_df[idx_in_mesh,]
    ML_cov_model_spat = ggplot(ML_cov_model_df) + geom_point(aes(x=x,y=y,color = space),size=5) + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    S_ef = S_ef + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    
    #grid.arrange(ML_cov_model_spat,S_ef)
    
    ##############
    ### SML cov
    ##############
    SML_cov_form <- y ~ 0 + b0 + cov_small  + cov_med + cov_large +
      f(spat, model=spde)   
    
    #### fit model
    SML_cov_model <- inla(  SML_cov_form,
                            family="poisson",
                            data=inla.stack.data(both),
                            control.compute=list(waic=T,cpo=T),
                            control.predictor=list( A=inla.stack.A(both), compute=TRUE, link=1),
                            verbose=FALSE,  #control.inla = ci,
                            num.threads = 1)
    
    SML_cov_model$summary.fixed
    
    SML_cov_model_df = data.frame(x=mesh$loc[,1],y=mesh$loc[,2],space = SML_cov_model$summary.random$spat$mean)
    SML_cov_model_df=SML_cov_model_df[idx_in_mesh,]
    SML_cov_model_spat = ggplot(SML_cov_model_df) + geom_point(aes(x=x,y=y,color = space),size=5) + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    resid = resid + geom_point(data=sample_data,aes(x=x,y=y),color="red",size=1)
    
    ##############
    ### only cov
    ##############
    only_cov_form <- y ~ 0 + b0 + cov_small  + cov_med + cov_large 
    
    #### fit model
    only_cov_model <- inla(  only_cov_form,
                             family="poisson",
                             data=inla.stack.data(both),
                             control.compute=list(waic=T,cpo=T),
                             control.predictor=list( A=inla.stack.A(both), compute=TRUE, link=1),
                             verbose=FALSE,  #control.inla = ci,
                             num.threads = 1)
    
    only_cov_model$summary.fixed
    
    
  }
  
  
  #####################################################
  ########  surfaces #############################
  ####################################################
  df.pred.spat = df.pred ; df.pred.spat[,3:5]=0
  
  fields = data.frame(abundance=abu_sp1,
                      x=round(coordinates(cov_small_surface)[,1],2),
                      y=round(coordinates(cov_small_surface)[,2],2),
                      ef_S=effect_small_cov,
                      ef_M=effect_med_cov,
                      ef_L=effect_large_cov,
                      ef_SL=effect_small_cov+effect_large_cov,
                      ef_SM=effect_med_cov+effect_small_cov,
                      ef_ML=effect_large_cov+effect_med_cov,
                      ef_SML = effect_small_cov+effect_large_cov+effect_med_cov,
                      
                      resid = scale(abu_sp1) - scale(effect_small_cov+effect_large_cov+effect_med_cov),
                      
                      M_0 =  no_cov_model$summary.random$spat$mean[idx_in_mesh],
                      M_S =  S_cov_model$summary.random$spat$mean[idx_in_mesh],
                      M_M =  M_cov_model$summary.random$spat$mean[idx_in_mesh],
                      M_L =  L_cov_model$summary.random$spat$mean[idx_in_mesh],
                      M_SM =  SM_cov_model$summary.random$spat$mean[idx_in_mesh],
                      M_SL =  SL_cov_model$summary.random$spat$mean[idx_in_mesh],
                      M_ML =  ML_cov_model$summary.random$spat$mean[idx_in_mesh],
                      M_SML =  SML_cov_model$summary.random$spat$mean[idx_in_mesh],
                      
                      M_0_gam =  predict(no_cov_gam,df.pred.spat,type = "response"),
                      M_S_gam =  predict(S_cov_gam,df.pred.spat,type = "response"),
                      M_M_gam =  predict(M_cov_gam,df.pred.spat,type = "response"),
                      M_L_gam =  predict(L_cov_gam,df.pred.spat,type = "response"),
                      M_SM_gam =  predict(SM_cov_gam,df.pred.spat,type = "response"),
                      M_SL_gam =  predict(SL_cov_gam,df.pred.spat,type = "response"),
                      M_ML_gam =  predict(ML_cov_gam,df.pred.spat,type = "response"),
                      M_SML_gam =  predict(SML_cov_gam,df.pred.spat,type = "response"),
                      
                      
                      M_only_cov_gam =  predict(only_cov_gam,df.pred,type = "response"),
                      M_only_cov_inla =  only_cov_model$summary.fitted.values[idx_pred,1][idx_in_mesh]
  )
  
  geo_points = SpatialPoints(coords=cbind(coordinates(cov_small_surface)[,1],coordinates(cov_small_surface)[,2]),
                             #data=fields,
                             proj4string = CRS(as.character(NA)))
  
  round(cor(fields),2)[4:11,12:19]
  
  
  m=SpatialPixelsDataFrame(geo_points, fields, tolerance = sqrt(.Machine$double.eps), 
                           proj4string = CRS(as.character(NA)))
  
  
  ###############
  ### SSIM #####
  ###############
  require(dplyr)
  require(SpatialPack)
  
  names(m)
  
  ra_abu = raster(m,layer=1)
  
  ra_cov_S = raster(m,layer=4)
  ra_cov_M = raster(m,layer=5)
  ra_cov_L = raster(m,layer=6)
  ra_cov_SL = raster(m,layer=7)
  ra_cov_SM = raster(m,layer=8)
  ra_cov_ML = raster(m,layer=9)
  ra_cov_SML = raster(m,layer=10)
  
  ra_resid = raster(m,layer=11)
  
  ra_fit_M_0 = raster(m,layer=12)
  ra_fit_M_S = raster(m,layer=13)
  ra_fit_M_M = raster(m,layer=14)
  ra_fit_M_L = raster(m,layer=15)
  ra_fit_M_SM = raster(m,layer=16)
  ra_fit_M_SL = raster(m,layer=17)
  ra_fit_M_ML = raster(m,layer=18)
  ra_fit_M_SML = raster(m,layer=19)
  
  ra_fit_M_0_gam = raster(m,layer=20)
  ra_fit_M_S_gam = raster(m,layer=21)
  ra_fit_M_M_gam = raster(m,layer=22)
  ra_fit_M_L_gam = raster(m,layer=23)
  ra_fit_M_SM_gam = raster(m,layer=24)
  ra_fit_M_SL_gam = raster(m,layer=25)
  ra_fit_M_ML_gam = raster(m,layer=26)
  ra_fit_M_SML_gam = raster(m,layer=27)
  
  ra_fit_M_only_cov = raster(m,layer=29)
  ra_fit_M_only_cov_gam = raster(m,layer=28)
  
  #### Calculate Similarity In Pattern (SIP)
  
  ####### M0 SIP comparisons ##########
  M_0_S_cov <- ssimMap(ra_cov_S,ra_fit_M_0)
  M_0_M_cov <- ssimMap(ra_cov_M,ra_fit_M_0)
  M_0_L_cov <- ssimMap(ra_cov_L,ra_fit_M_0)
  M_0_SM_cov <- ssimMap(ra_cov_SM,ra_fit_M_0)
  M_0_SL_cov <- ssimMap(ra_cov_SL,ra_fit_M_0)
  M_0_ML_cov <- ssimMap(ra_cov_ML,ra_fit_M_0)
  M_0_SML_cov <- ssimMap(ra_cov_SML,ra_fit_M_0)
  M_0_resid <- ssimMap(ra_resid,ra_fit_M_0)
  
  
  M_0_SIP = c(
    median(raster::values(M_0_S_cov$SIP),na.rm=T),
    median(raster::values(M_0_M_cov$SIP),na.rm=T),
    median(raster::values(M_0_L_cov$SIP),na.rm=T),
    median(raster::values(M_0_SM_cov$SIP),na.rm=T),
    median(raster::values(M_0_SL_cov$SIP),na.rm=T),
    median(raster::values(M_0_ML_cov$SIP),na.rm=T),
    median(raster::values(M_0_SML_cov$SIP),na.rm=T),
    median(raster::values(M_0_resid$SIP),na.rm=T)
  )
  
  
  M_0_S_cov_gam <- ssimMap(ra_cov_S,ra_fit_M_0_gam)
  M_0_M_cov_gam <- ssimMap(ra_cov_M,ra_fit_M_0_gam)
  M_0_L_cov_gam <- ssimMap(ra_cov_L,ra_fit_M_0_gam)
  M_0_SM_cov_gam <- ssimMap(ra_cov_SM,ra_fit_M_0_gam)
  M_0_SL_cov_gam <- ssimMap(ra_cov_SL,ra_fit_M_0_gam)
  M_0_ML_cov_gam <- ssimMap(ra_cov_ML,ra_fit_M_0_gam)
  M_0_SML_cov_gam <- ssimMap(ra_cov_SML,ra_fit_M_0_gam)
  M_0_resid_gam <- ssimMap(ra_resid,ra_fit_M_0_gam)
  
  
  M_0_SIP_gam = c(
    median(raster::values(M_0_S_cov_gam$SIP),na.rm=T),
    median(raster::values(M_0_M_cov_gam$SIP),na.rm=T),
    median(raster::values(M_0_L_cov_gam$SIP),na.rm=T),
    median(raster::values(M_0_SM_cov_gam$SIP),na.rm=T),
    median(raster::values(M_0_SL_cov_gam$SIP),na.rm=T),
    median(raster::values(M_0_ML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_0_SML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_0_resid_gam$SIP),na.rm=T)
  )
  
  ####### M_S SIP comparisons ##########
  M_S_S_cov <- ssimMap(ra_cov_S,ra_fit_M_S)
  M_S_M_cov <- ssimMap(ra_cov_M,ra_fit_M_S)
  M_S_L_cov <- ssimMap(ra_cov_L,ra_fit_M_S)
  M_S_SM_cov <- ssimMap(ra_cov_SM,ra_fit_M_S)
  M_S_SL_cov <- ssimMap(ra_cov_SL,ra_fit_M_S)
  M_S_ML_cov <- ssimMap(ra_cov_ML,ra_fit_M_S)
  M_S_SML_cov <- ssimMap(ra_cov_SML,ra_fit_M_S)
  M_S_resid <- ssimMap(ra_resid,ra_fit_M_S)
  
  
  M_S_SIP = c(
    median(raster::values(M_S_S_cov$SIP),na.rm=T),
    median(raster::values(M_S_M_cov$SIP),na.rm=T),
    median(raster::values(M_S_L_cov$SIP),na.rm=T),
    median(raster::values(M_S_SM_cov$SIP),na.rm=T),
    median(raster::values(M_S_SL_cov$SIP),na.rm=T),
    median(raster::values(M_S_ML_cov$SIP),na.rm=T),
    median(raster::values(M_S_SML_cov$SIP),na.rm=T),
    median(raster::values(M_S_resid$SIP),na.rm=T)
  )
  
  M_S_S_cov_gam <- ssimMap(ra_cov_S,ra_fit_M_S_gam)
  M_S_M_cov_gam <- ssimMap(ra_cov_M,ra_fit_M_S_gam)
  M_S_L_cov_gam <- ssimMap(ra_cov_L,ra_fit_M_S_gam)
  M_S_SM_cov_gam <- ssimMap(ra_cov_SM,ra_fit_M_S_gam)
  M_S_SL_cov_gam <- ssimMap(ra_cov_SL,ra_fit_M_S_gam)
  M_S_ML_cov_gam <- ssimMap(ra_cov_ML,ra_fit_M_S_gam)
  M_S_SML_cov_gam <- ssimMap(ra_cov_SML,ra_fit_M_S_gam)
  M_S_resid_gam <- ssimMap(ra_resid,ra_fit_M_S_gam)
  
  
  M_S_SIP_gam = c(
    median(raster::values(M_S_S_cov_gam$SIP),na.rm=T),
    median(raster::values(M_S_M_cov_gam$SIP),na.rm=T),
    median(raster::values(M_S_L_cov_gam$SIP),na.rm=T),
    median(raster::values(M_S_SM_cov_gam$SIP),na.rm=T),
    median(raster::values(M_S_SL_cov_gam$SIP),na.rm=T),
    median(raster::values(M_S_ML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_S_SML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_S_resid_gam$SIP),na.rm=T)
  )
  
  ####### M_M SIP comparisons ##########
  M_M_S_cov <- ssimMap(ra_cov_S,ra_fit_M_M)
  M_M_M_cov <- ssimMap(ra_cov_M,ra_fit_M_M)
  M_M_L_cov <- ssimMap(ra_cov_L,ra_fit_M_M)
  M_M_SM_cov <- ssimMap(ra_cov_SM,ra_fit_M_M)
  M_M_SL_cov <- ssimMap(ra_cov_SL,ra_fit_M_M)
  M_M_ML_cov <- ssimMap(ra_cov_ML,ra_fit_M_M)
  M_M_SML_cov <- ssimMap(ra_cov_SML,ra_fit_M_M)
  M_M_resid <- ssimMap(ra_resid,ra_fit_M_M)
  
  
  M_M_SIP = c(
    median(raster::values(M_M_S_cov$SIP),na.rm=T),
    median(raster::values(M_M_M_cov$SIP),na.rm=T),
    median(raster::values(M_M_L_cov$SIP),na.rm=T),
    median(raster::values(M_M_SM_cov$SIP),na.rm=T),
    median(raster::values(M_M_SL_cov$SIP),na.rm=T),
    median(raster::values(M_M_ML_cov$SIP),na.rm=T),
    median(raster::values(M_M_SML_cov$SIP),na.rm=T),
    median(raster::values(M_M_resid$SIP),na.rm=T)
  )
  
  M_M_S_cov_gam <- ssimMap(ra_cov_S,ra_fit_M_M_gam)
  M_M_M_cov_gam <- ssimMap(ra_cov_M,ra_fit_M_M_gam)
  M_M_L_cov_gam <- ssimMap(ra_cov_L,ra_fit_M_M_gam)
  M_M_SM_cov_gam <- ssimMap(ra_cov_SM,ra_fit_M_M_gam)
  M_M_SL_cov_gam <- ssimMap(ra_cov_SL,ra_fit_M_M_gam)
  M_M_ML_cov_gam <- ssimMap(ra_cov_ML,ra_fit_M_M_gam)
  M_M_SML_cov_gam <- ssimMap(ra_cov_SML,ra_fit_M_M_gam)
  M_M_resid_gam <- ssimMap(ra_resid,ra_fit_M_M_gam)
  
  
  M_M_SIP_gam = c(
    median(raster::values(M_M_S_cov_gam$SIP),na.rm=T),
    median(raster::values(M_M_M_cov_gam$SIP),na.rm=T),
    median(raster::values(M_M_L_cov_gam$SIP),na.rm=T),
    median(raster::values(M_M_SM_cov_gam$SIP),na.rm=T),
    median(raster::values(M_M_SL_cov_gam$SIP),na.rm=T),
    median(raster::values(M_M_ML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_M_SML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_M_resid_gam$SIP),na.rm=T)
  )
  
  ####### ML SIP comparisons ##########
  M_L_S_cov <- ssimMap(ra_cov_S,ra_fit_M_L)
  M_L_M_cov <- ssimMap(ra_cov_M,ra_fit_M_L)
  M_L_L_cov <- ssimMap(ra_cov_L,ra_fit_M_L)
  M_L_SM_cov <- ssimMap(ra_cov_SM,ra_fit_M_L)
  M_L_SL_cov <- ssimMap(ra_cov_SL,ra_fit_M_L)
  M_L_ML_cov <- ssimMap(ra_cov_ML,ra_fit_M_L)
  M_L_SML_cov <- ssimMap(ra_cov_SML,ra_fit_M_L)
  M_L_resid <- ssimMap(ra_resid,ra_fit_M_L)
  
  
  M_L_SIP = c(
    median(raster::values(M_L_S_cov$SIP),na.rm=T),
    median(raster::values(M_L_M_cov$SIP),na.rm=T),
    median(raster::values(M_L_L_cov$SIP),na.rm=T),
    median(raster::values(M_L_SM_cov$SIP),na.rm=T),
    median(raster::values(M_L_SL_cov$SIP),na.rm=T),
    median(raster::values(M_L_ML_cov$SIP),na.rm=T),
    median(raster::values(M_L_SML_cov$SIP),na.rm=T),
    median(raster::values(M_L_resid$SIP),na.rm=T)
  )
  
  M_L_S_cov_gam <- ssimMap(ra_cov_S,ra_fit_M_L_gam)
  M_L_M_cov_gam <- ssimMap(ra_cov_M,ra_fit_M_L_gam)
  M_L_L_cov_gam <- ssimMap(ra_cov_L,ra_fit_M_L_gam)
  M_L_SM_cov_gam <- ssimMap(ra_cov_SM,ra_fit_M_L_gam)
  M_L_SL_cov_gam <- ssimMap(ra_cov_SL,ra_fit_M_L_gam)
  M_L_ML_cov_gam <- ssimMap(ra_cov_ML,ra_fit_M_L_gam)
  M_L_SML_cov_gam <- ssimMap(ra_cov_SML,ra_fit_M_L_gam)
  M_L_resid_gam <- ssimMap(ra_resid,ra_fit_M_L_gam)
  
  
  M_L_SIP_gam = c(
    median(raster::values(M_L_S_cov_gam$SIP),na.rm=T),
    median(raster::values(M_L_M_cov_gam$SIP),na.rm=T),
    median(raster::values(M_L_L_cov_gam$SIP),na.rm=T),
    median(raster::values(M_L_SM_cov_gam$SIP),na.rm=T),
    median(raster::values(M_L_SL_cov_gam$SIP),na.rm=T),
    median(raster::values(M_L_ML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_L_SML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_L_resid_gam$SIP),na.rm=T)
  )
  
  ####### MSM SIP comparisons ##########
  M_SM_S_cov <- ssimMap(ra_cov_S,ra_fit_M_SM)
  M_SM_M_cov <- ssimMap(ra_cov_M,ra_fit_M_SM)
  M_SM_L_cov <- ssimMap(ra_cov_L,ra_fit_M_SM)
  M_SM_SM_cov <- ssimMap(ra_cov_SM,ra_fit_M_SM)
  M_SM_SL_cov <- ssimMap(ra_cov_SL,ra_fit_M_SM)
  M_SM_ML_cov <- ssimMap(ra_cov_ML,ra_fit_M_SM)
  M_SM_SML_cov <- ssimMap(ra_cov_SML,ra_fit_M_SM)
  M_SM_resid <- ssimMap(ra_resid,ra_fit_M_SM)
  
  
  M_SM_SIP = c(
    median(raster::values(M_SM_S_cov$SIP),na.rm=T),
    median(raster::values(M_SM_M_cov$SIP),na.rm=T),
    median(raster::values(M_SM_L_cov$SIP),na.rm=T),
    median(raster::values(M_SM_SM_cov$SIP),na.rm=T),
    median(raster::values(M_SM_SL_cov$SIP),na.rm=T),
    median(raster::values(M_SM_ML_cov$SIP),na.rm=T),
    median(raster::values(M_SM_SML_cov$SIP),na.rm=T),
    median(raster::values(M_SM_resid$SIP),na.rm=T)
  )
  
  M_SM_S_cov_gam <- ssimMap(ra_cov_S,ra_fit_M_SM_gam)
  M_SM_M_cov_gam <- ssimMap(ra_cov_M,ra_fit_M_SM_gam)
  M_SM_L_cov_gam <- ssimMap(ra_cov_L,ra_fit_M_SM_gam)
  M_SM_SM_cov_gam <- ssimMap(ra_cov_SM,ra_fit_M_SM_gam)
  M_SM_SL_cov_gam <- ssimMap(ra_cov_SL,ra_fit_M_SM_gam)
  M_SM_ML_cov_gam <- ssimMap(ra_cov_ML,ra_fit_M_SM_gam)
  M_SM_SML_cov_gam <- ssimMap(ra_cov_SML,ra_fit_M_SM_gam)
  M_SM_resid_gam <- ssimMap(ra_resid,ra_fit_M_SM_gam)
  
  
  M_SM_SIP_gam = c(
    median(raster::values(M_SM_S_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SM_M_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SM_L_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SM_SM_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SM_SL_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SM_ML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SM_SML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SM_resid_gam$SIP),na.rm=T)
  )
  
  ####### MSL SIP comparisons ##########
  M_SL_S_cov <- ssimMap(ra_cov_S,ra_fit_M_SL)
  M_SL_M_cov <- ssimMap(ra_cov_M,ra_fit_M_SL)
  M_SL_L_cov <- ssimMap(ra_cov_L,ra_fit_M_SL)
  M_SL_SM_cov <- ssimMap(ra_cov_SM,ra_fit_M_SL)
  M_SL_SL_cov <- ssimMap(ra_cov_SL,ra_fit_M_SL)
  M_SL_ML_cov <- ssimMap(ra_cov_ML,ra_fit_M_SL)
  M_SL_SML_cov <- ssimMap(ra_cov_SML,ra_fit_M_SL)
  M_SL_resid <- ssimMap(ra_resid,ra_fit_M_SL)
  
  
  M_SL_SIP = c(
    median(raster::values(M_SL_S_cov$SIP),na.rm=T),
    median(raster::values(M_SL_M_cov$SIP),na.rm=T),
    median(raster::values(M_SL_L_cov$SIP),na.rm=T),
    median(raster::values(M_SL_SM_cov$SIP),na.rm=T),
    median(raster::values(M_SL_SL_cov$SIP),na.rm=T),
    median(raster::values(M_SL_ML_cov$SIP),na.rm=T),
    median(raster::values(M_SL_SML_cov$SIP),na.rm=T),
    median(raster::values(M_SL_resid$SIP),na.rm=T)
  )
  
  M_SL_S_cov_gam <- ssimMap(ra_cov_S,ra_fit_M_SL_gam)
  M_SL_M_cov_gam <- ssimMap(ra_cov_M,ra_fit_M_SL_gam)
  M_SL_L_cov_gam <- ssimMap(ra_cov_L,ra_fit_M_SL_gam)
  M_SL_SM_cov_gam <- ssimMap(ra_cov_SM,ra_fit_M_SL_gam)
  M_SL_SL_cov_gam <- ssimMap(ra_cov_SL,ra_fit_M_SL_gam)
  M_SL_ML_cov_gam <- ssimMap(ra_cov_ML,ra_fit_M_SL_gam)
  M_SL_SML_cov_gam <- ssimMap(ra_cov_SML,ra_fit_M_SL_gam)
  M_SL_resid_gam <- ssimMap(ra_resid,ra_fit_M_SL_gam)
  
  
  M_SL_SIP_gam = c(
    median(raster::values(M_SL_S_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SL_M_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SL_L_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SL_SM_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SL_SL_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SL_ML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SL_SML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SL_resid_gam$SIP),na.rm=T)
  )
  
  ####### MML SIP comparisons ##########
  M_ML_S_cov <- ssimMap(ra_cov_S,ra_fit_M_ML)
  M_ML_M_cov <- ssimMap(ra_cov_M,ra_fit_M_ML)
  M_ML_L_cov <- ssimMap(ra_cov_L,ra_fit_M_ML)
  M_ML_SM_cov <- ssimMap(ra_cov_SM,ra_fit_M_ML)
  M_ML_SL_cov <- ssimMap(ra_cov_SL,ra_fit_M_ML)
  M_ML_ML_cov <- ssimMap(ra_cov_ML,ra_fit_M_ML)
  M_ML_SML_cov <- ssimMap(ra_cov_SML,ra_fit_M_ML)
  M_ML_resid <- ssimMap(ra_resid,ra_fit_M_ML)
  
  
  M_ML_SIP = c(
    median(raster::values(M_ML_S_cov$SIP),na.rm=T),
    median(raster::values(M_ML_M_cov$SIP),na.rm=T),
    median(raster::values(M_ML_L_cov$SIP),na.rm=T),
    median(raster::values(M_ML_SM_cov$SIP),na.rm=T),
    median(raster::values(M_ML_SL_cov$SIP),na.rm=T),
    median(raster::values(M_ML_ML_cov$SIP),na.rm=T),
    median(raster::values(M_ML_SML_cov$SIP),na.rm=T),
    median(raster::values(M_ML_resid$SIP),na.rm=T)
  )
  
  M_ML_S_cov_gam <- ssimMap(ra_cov_S,ra_fit_M_ML_gam)
  M_ML_M_cov_gam <- ssimMap(ra_cov_M,ra_fit_M_ML_gam)
  M_ML_L_cov_gam <- ssimMap(ra_cov_L,ra_fit_M_ML_gam)
  M_ML_SM_cov_gam <- ssimMap(ra_cov_SM,ra_fit_M_ML_gam)
  M_ML_SL_cov_gam <- ssimMap(ra_cov_SL,ra_fit_M_ML_gam)
  M_ML_ML_cov_gam <- ssimMap(ra_cov_ML,ra_fit_M_ML_gam)
  M_ML_SML_cov_gam <- ssimMap(ra_cov_SML,ra_fit_M_ML_gam)
  M_ML_resid_gam <- ssimMap(ra_resid,ra_fit_M_ML_gam)
  
  
  M_ML_SIP_gam = c(
    median(raster::values(M_ML_S_cov_gam$SIP),na.rm=T),
    median(raster::values(M_ML_M_cov_gam$SIP),na.rm=T),
    median(raster::values(M_ML_L_cov_gam$SIP),na.rm=T),
    median(raster::values(M_ML_SM_cov_gam$SIP),na.rm=T),
    median(raster::values(M_ML_SL_cov_gam$SIP),na.rm=T),
    median(raster::values(M_ML_ML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_ML_SML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_ML_resid_gam$SIP),na.rm=T)
  )
  
  ####### MSML SIP comparisons ##########
  M_SML_S_cov <- ssimMap(ra_cov_S,ra_fit_M_SML)
  M_SML_M_cov <- ssimMap(ra_cov_M,ra_fit_M_SML)
  M_SML_L_cov <- ssimMap(ra_cov_L,ra_fit_M_SML)
  M_SML_SM_cov <- ssimMap(ra_cov_SM,ra_fit_M_SML)
  M_SML_SL_cov <- ssimMap(ra_cov_SL,ra_fit_M_SML)
  M_SML_ML_cov <- ssimMap(ra_cov_ML,ra_fit_M_SML)
  M_SML_SML_cov <- ssimMap(ra_cov_SML,ra_fit_M_SML)
  M_SML_resid <- ssimMap(ra_resid,ra_fit_M_SML)
  
  
  M_SML_SIP = c(
    median(raster::values(M_SML_S_cov$SIP),na.rm=T),
    median(raster::values(M_SML_M_cov$SIP),na.rm=T),
    median(raster::values(M_SML_L_cov$SIP),na.rm=T),
    median(raster::values(M_SML_SM_cov$SIP),na.rm=T),
    median(raster::values(M_SML_SL_cov$SIP),na.rm=T),
    median(raster::values(M_SML_ML_cov$SIP),na.rm=T),
    median(raster::values(M_SML_SML_cov$SIP),na.rm=T),
    median(raster::values(M_SML_resid$SIP),na.rm=T)
  )
  
  M_SML_S_cov_gam <- ssimMap(ra_cov_S,ra_fit_M_SML_gam)
  M_SML_M_cov_gam <- ssimMap(ra_cov_M,ra_fit_M_SML_gam)
  M_SML_L_cov_gam <- ssimMap(ra_cov_L,ra_fit_M_SML_gam)
  M_SML_SM_cov_gam <- ssimMap(ra_cov_SM,ra_fit_M_SML_gam)
  M_SML_SL_cov_gam <- ssimMap(ra_cov_SL,ra_fit_M_SML_gam)
  M_SML_ML_cov_gam <- ssimMap(ra_cov_ML,ra_fit_M_SML_gam)
  M_SML_SML_cov_gam <- ssimMap(ra_cov_SML,ra_fit_M_SML_gam)
  M_SML_resid_gam <- ssimMap(ra_resid,ra_fit_M_SML_gam)
  
  
  M_SML_SIP_gam = c(
    median(raster::values(M_SML_S_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SML_M_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SML_L_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SML_SM_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SML_SL_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SML_ML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SML_SML_cov_gam$SIP),na.rm=T),
    median(raster::values(M_SML_resid_gam$SIP),na.rm=T)
  )
  
  SIP_matrix = round(rbind(M_0_SIP,M_S_SIP,M_M_SIP,M_L_SIP,M_SM_SIP,M_SL_SIP,M_ML_SIP,M_SML_SIP),2)
  delta_matrix = round(rbind(abs(max(M_0_SIP)-M_0_SIP),abs(max(M_S_SIP)-M_S_SIP),
                             abs(max(M_M_SIP)-M_M_SIP),abs(max(M_L_SIP)-M_L_SIP),
                             abs(max(M_SM_SIP)-M_SM_SIP),abs(max(M_SL_SIP)-M_SL_SIP),
                             abs(max(M_ML_SIP)-M_ML_SIP),abs(max(M_SML_SIP)-M_SML_SIP)),2)
  
  
  delta_M_0_SIP[sim,]= abs(max(M_0_SIP)-M_0_SIP)
  delta_M_0_SIP_gam[sim,]= abs(max(M_0_SIP_gam)-M_0_SIP_gam)
  
  delta_M_S_SIP[sim,]= abs(max(M_S_SIP)-M_S_SIP)
  delta_M_S_SIP_gam[sim,]= abs(max(M_S_SIP_gam)-M_S_SIP_gam)
  
  delta_M_M_SIP[sim,]= abs(max(M_M_SIP)-M_M_SIP)
  delta_M_M_SIP_gam[sim,]= abs(max(M_M_SIP_gam)-M_M_SIP_gam)
  
  delta_M_L_SIP[sim,]= abs(max(M_L_SIP)-M_L_SIP)
  delta_M_L_SIP_gam[sim,]= abs(max(M_L_SIP_gam)-M_L_SIP_gam)
  
  delta_M_SM_SIP[sim,]= abs(max(M_SM_SIP)-M_SM_SIP)
  delta_M_SM_SIP_gam[sim,]= abs(max(M_SM_SIP_gam)-M_SM_SIP_gam)
  
  delta_M_SL_SIP[sim,]= abs(max(M_SL_SIP)-M_SL_SIP)
  delta_M_SL_SIP_gam[sim,]= abs(max(M_SL_SIP_gam)-M_SL_SIP_gam)
  
  delta_M_ML_SIP[sim,]= abs(max(M_ML_SIP)-M_ML_SIP)
  delta_M_ML_SIP_gam[sim,]= abs(max(M_ML_SIP_gam)-M_ML_SIP_gam)
  
  delta_M_SML_SIP[sim,]= abs(max(M_SML_SIP)-M_SML_SIP)
  delta_M_SML_SIP_gam[sim,]= abs(max(M_SML_SIP_gam)-M_SML_SIP_gam)
  
  
  
  print(paste("sim =",sim, round(Sys.time()-time0,2)))
  
}



colnames(delta_M_0_SIP) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_S_SIP) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_M_SIP) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_L_SIP) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_SM_SIP) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_SL_SIP) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_ML_SIP) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_SML_SIP) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  

colnames(delta_M_0_SIP_gam) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_S_SIP_gam) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_M_SIP_gam) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_L_SIP_gam) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_SM_SIP_gam) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_SL_SIP_gam) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_ML_SIP_gam) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  
colnames(delta_M_SML_SIP_gam) = c("Cov_S","Cov_M","Cov_L","Cov_SM","Cov_SL","Cov_ML","Cov_SML","Resid")  



####################################
########## INLA results #########
#####################################

SIP_results = list(delta_M_0_SIP=delta_M_0_SIP,delta_M_S_SIP=delta_M_S_SIP,delta_M_M_SIP=delta_M_M_SIP,delta_M_L_SIP=delta_M_L_SIP,
                   delta_M_SM_SIP=delta_M_SM_SIP,delta_M_SL_SIP=delta_M_SL_SIP,delta_M_ML_SIP=delta_M_ML_SIP,
                   delta_M_SML_SIP=delta_M_SML_SIP)


SIP_results_gam = list(delta_M_0_SIP=delta_M_0_SIP_gam,delta_M_S_SIP=delta_M_S_SIP_gam,delta_M_M_SIP=delta_M_M_SIP_gam,delta_M_L_SIP=delta_M_L_SIP_gam,
                       delta_M_SM_SIP=delta_M_SM_SIP_gam,delta_M_SL_SIP=delta_M_SL_SIP_gam,delta_M_ML_SIP=delta_M_ML_SIP_gam,
                       delta_M_SML_SIP=delta_M_SML_SIP_gam)




####################################
########## MGCV results #########
#####################################
SIP_results_mean = SIP_results_median = SIP_results_sd = matrix(NA,ncol=ncol(delta_mae)-1,nrow=length(SIP_results_gam))
for(i in 1:length(SIP_results_gam)){
  SIP_results_mean[i,] = round(apply(SIP_results_gam[[i]],2,mean),2)
  SIP_results_median[i,] = round(apply(SIP_results_gam[[i]],2,median),2)
  SIP_results_sd[i,] = round(apply(SIP_results_gam[[i]],2,sd),2)
}
row.names(SIP_results_mean) = row.names(SIP_results_median) = row.names(SIP_results_sd) = names(SIP_results)
colnames(SIP_results_mean) = colnames(SIP_results_median) = colnames(SIP_results_sd) = colnames(delta_M_0_SIP)

xtable(SIP_results_mean)
xtable(SIP_results_sd)


score_results_mean = score_results_median = score_results_sd = matrix(NA,ncol=ncol(delta_mae),nrow=length(score_results_gam))
for(i in 1:length(score_results_gam)){
  score_results_mean[i,] = apply(score_results_gam[[i]],2,mean)
  score_results_median[i,] = apply(score_results_gam[[i]],2,median)
  score_results_sd[i,] = apply(score_results_gam[[i]],2,sd)
}
row.names(score_results_mean) = row.names(score_results_median) = row.names(score_results_sd) = names(score_results_gam)
colnames(score_results_mean) = colnames(score_results_median) = colnames(score_results_sd) = colnames(delta_mae)



