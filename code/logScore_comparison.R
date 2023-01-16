
#####   compute log score comparison for given training and test data   #####

compScore=function(locs,data.train,data.test,methods=1:4,NNarray.max,
                   cor.ord=FALSE,m.max=30,holdout=FALSE){

  ## general settings and parameters
  n=nrow(data.train)
  N=ncol(data.train)
  N.test=ncol(data.test)
  log.scores=array(dim=c(length(methods),N.test))
  if(holdout) ho.scores=array(dim=c(length(methods),N.test))
  n.fix=round(n*.5)
  k=1  # method counter

  ## ordering and conditioning
  if(missing(NNarray.max)){
    if(cor.ord) {
      dists=rdist(locs)
      cor.est=cor(t(data.train))*exp(-dists/max(dists))
      ord=order_maxmin_correlation(cor.est)
      NNarray.max=find_nn(cor.est[ord,ord],m.max)
      scales=computeScales_cor(cor.est[ord,ord],NNarray.max)
    } else {
      ord=GPvecchia::order_maxmin_exact(locs)
      NNarray.max=GpGp::find_ordered_nn(locs[ord,],m.max)[,-1]
      scales=computeScales(locs[ord,],NNarray.max)
    }
    dat.train=data.train[ord,]
    dat.test=data.test[ord,]
  } else {
    dat.train=data.train
    dat.test=data.test
    scales=computeScales(locs,NNarray.max)
  }

  
  ## serial or parallel depending on N
  parallel = (N>=100 & n>=500)  ######################
  
  ## 1: linear (w/ UQ)
  if(1 %in% methods){ 
    
    print('linear')
    
    # train model
    fit=optimFitMap(dat.train,NNarray.max,scales=scales,
                    linear=TRUE,parallel=parallel)

    # compute log scores
    for(j in 1:N.test) log.scores[k,j]=condSamp(fit,mode='score',obs=dat.test[,j])
    
    # compute other scores if requested
    if(holdout) ho.scores[k,]=compHoldout(fit,dat.test,mode='score',n.fix)
    
    k=k+1
    
  }
  
  
  ## 2: linear (w/o UQ)
  if(2 %in% methods){ 
    
    print('linear (w/o UQ)')
    
    # compute log scores
    for(j in 1:N.test) log.scores[k,j]=condSamp(fit,mode='scorepm',
                                                obs=dat.test[,j])
    
    # compute other scores if requested
    if(holdout) ho.scores[k,]=compHoldout(fit,dat.test,mode='scorepm',n.fix)
    
    k=k+1
    
  }
  
  
  ## 3: nonlinear (w/ UQ)
  if(3 %in% methods){ 
    
    print('nonlinear')
    
    # train model
    fit=optimFitMap(dat.train,NNarray.max,scales=scales,
                    linear=FALSE,parallel=parallel)
    
    # compute log scores
    for(j in 1:N.test) log.scores[k,j]=condSamp(fit,mode='score',obs=dat.test[,j])
    
    # compute other scores if requested
    if(holdout) ho.scores[k,]=compHoldout(fit,dat.test,mode='score',n.fix)
    
    k=k+1
    
  }

      
  ## 4: nonlinear (w/o UQ)
  if(4 %in% methods){ 
    
    print('nonlinear (w/o UQ)')
    
    # compute log scores
    for(j in 1:N.test) log.scores[k,j]=condSamp(fit,mode='scorepm',
                                                obs=dat.test[,j])
    
    # compute other scores if requested
    if(holdout) ho.scores[k,]=compHoldout(fit,dat.test,mode='scorepm',n.fix)
    
    k=k+1
    
  }
    
    
  ## 5: DPM + nonlinear
  if(5 %in% methods){ 
    
    print('DPM nonlinear')
    
    # train model
    fit=fitDPM(dat.train,NNarray.max,L=500,burnin=200,thin.by=10,
               scales=scales,parallel=TRUE)
    
    # compute log scores
    for(j in 1:N.test) log.scores[k,j]=predDPM(fit,obs=dat.test[,j],
                                               parallel=TRUE,mode='score')
    
    if(holdout){
      for(j in 1:N.test) ho.scores[k,j]=predDPMSerial(fit,mode='score',
                                obs=dat.test[,j],x.fixed=dat.test[1:n.fix,j])
    }
    
    k=k+1
    
  }
  
  
  ## return values
  if(holdout){
    hoi=ord[(n.fix+1):n]
    return(list(log.scores=log.scores,ho.scores=ho.scores,hoi=hoi)) 
  } else return(log.scores)
  
}


## function to compute LS for hold-out set
compHoldout=function(fit,dat.test,mode='score',n.fix){
  
  log.score=numeric(length=N.test)
  for(j in 1:N.test)
    log.score[j]=condSamp(fit,mode=mode,obs=dat.test[,j],
                                           x.fixed=dat.test[1:n.fix,j])
  
  return(log.score)
}




#################   other methods with which to compare    ##############

library(mvtnorm)


#### (gaussian) LS for hold-out set
ls.ho.gaussian=function(sigma,data.test,holdout){
  temp=solve(sigma[-holdout,-holdout],sigma[-holdout,holdout])
  postcov=sigma[holdout,holdout]-t(temp)%*%sigma[-holdout,holdout]
  postmean=t(temp)%*%data.test[-holdout,]
  ls=numeric(length=ncol(data.test))
  for(j in 1:length(ls)) ls[j]=mvtnorm::dmvnorm(data.test[holdout,j],
                                                postmean[,j],postcov,log=TRUE)
  mean(ls)
}

#### functions for GP with matern cov
matern_nloglik=function(params,dat,dists){
  cov.mat=exp(params[1])*fields::Matern(dists,exp(params[2]),nu=exp(params[3]))
  -mean(mvtnorm::dmvnorm(t(dat),rep(0,n),sigma=cov.mat,log=TRUE))
}
matern_ls_ho=function(params,dat,dists,holdout){
  cov.mat=exp(params[1])*fields::Matern(dists,exp(params[2]),nu=exp(params[3]))
  ls.ho.gaussian(cov.mat,dat,holdout)
}
matern_scores=function(locs,data.train,data.test,holdout=NULL){
  dists=rdist(locs)
  opt=optim(log(c(1,max(dists)*.1,1)),matern_nloglik,dat=data.train,dists=dists,
            control=list(trace=0,maxit=100,reltol=1e-3))
  ls=-matern_nloglik(opt$par,data.test,dists)
  if(is.null(holdout)) return(ls) else {
    ho.scores=matern_ls_ho(opt$par,data.test,dists,holdout)
    return(c(ls,mean(ho.scores)))
  }
}


####  locally stationary matern

# devtools::install_github("ashtonwiens/nonstationary")

## create local windows
local_windows=function(locs,ws=4){
  # ws is window size (how many in each direction)
  
  #get values from the grid
  lon=locs[,1]; lat=locs[,2]
  x <- unique(lon)
  y <- unique(lat)
  sizx <- diff(unique(lon))[1]
  sizy <- diff(unique(lat))[1]
  lon_lat <- cbind(lon, lat)
  grd2d_ord <- lon_lat
  nx <- length(x)
  ny <- length(y)
  
  #get neighbors in reduced window size
  ws=2
  windows_w <- matrix(NA, nr=nrow(grd2d_ord), nc=(ws*2+1)^2)
  for(i in 1:nrow(grd2d_ord)){
    temp <- grd2d_ord[i,]
    temp2 <- which(grd2d_ord[,1] <= (temp[[1]] + ws*sizx) & 
                     grd2d_ord[,1] >= (temp[[1]] - ws*sizx) &
                     grd2d_ord[,2] >= (temp[[2]] - ws*sizy) & 
                     grd2d_ord[,2] <= (temp[[2]] + ws*sizy)) 
    # store neighbors in window
    windows_w[i, 1:length(temp2)] <- temp2 
  }
  
  return(list(grd2d=grd2d_ord,windows=windows_w,nx=nx,ny=ny,x=x,y=y))
  
}


#### compute local matern scores

local_scores=function(grdwin,data.train,data.test,holdout=NULL){
  
  nx=grdwin$nx; ny=grdwin$ny
  
  ## fit local parameters for each window
  params <- matrix(NA, nc=5, nr=nx*ny)
  for(i in 1:(nx*ny)){
    sink('output/trash')
    huh <- nonstationary::fit.Matern.aniso(grdwin$grd2d[na.omit(grdwin$windows[i,]),], 
                            (data.train[na.omit(grdwin$windows[i,]),]))
    sink()
    params[i,] <- as.vector(huh$par)
  }
  
  ## parameter renormalization
  Lx <- matrix(sqrt(exp(params[,1])), nr = nx, nc = ny)
  Ly <- matrix(sqrt(exp(params[,4])),  nr = nx, nc = ny)
  sig <- matrix(sqrt(exp(params[,2])), nr = nx, nc = ny)
  tau2 <- matrix(exp(params[,3]),  nr = nx, nc = ny)
  angl <- matrix(params[,5],  nr = nx, nc = ny)
  angl <- (pi/2)*nonstationary::ilogit(angl) - (pi/4)
  
  ## combine all the local parameters
  testc <- nonstationary::Anisotropic.Matern.Nonstationary.Cov(grdwin$x, grdwin$y, 
                                Lx, Ly, angl, nu=1, sigma=sig) + diag(c(tau2))

  ## log likelihood calculation
  ls=mean(mvtnorm::dmvnorm(t(data.test), rep(0,n), sigma=testc, log=TRUE))
  
  if(is.null(holdout)) return(ls) else {
    ho.scores=ls.ho.gaussian(testc,data.test,holdout)
    return(c(ls,mean(ho.scores)))
  }

}
