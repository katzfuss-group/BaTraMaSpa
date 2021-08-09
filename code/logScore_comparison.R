
#####   compute log score comparison for given training and test data   #####

compScore=function(locs,data.train,data.test,methods=1:4,NNarray.max,
                   cor.ord=FALSE,m.max=30){

  ## general settings and parameters
  n=nrow(data.train)
  N=ncol(data.train)
  N.test=ncol(data.test)
  log.scores=array(dim=c(length(methods),N.test))
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
  parallel = (N>=100 & n>=500)
  
  ## 1: linear (w/ UQ)
  if(1 %in% methods){ 
    
    print('linear')
    
    # train model
    fit=optimFitMap(dat.train,NNarray.max,scales=scales,
                    linear=TRUE,parallel=parallel)

    # compute log scores
    for(j in 1:N.test) log.scores[k,j]=condSamp(fit,mode='score',obs=dat.test[,j])
    k=k+1
    
  }
  
  
  ## 2: linear (w/o UQ)
  if(2 %in% methods){ 
    
    print('linear (w/o UQ)')
    
    # compute log scores
    for(j in 1:N.test) log.scores[k,j]=condSamp(fit,mode='scorepm',
                                                obs=dat.test[,j])
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
    k=k+1
    
  }

      
  ## 4: nonlinear (w/o UQ)
  if(4 %in% methods){ 
    
    print('nonlinear (w/o UQ)')
    
    # compute log scores
    for(j in 1:N.test) log.scores[k,j]=condSamp(fit,mode='scorepm',
                                                obs=dat.test[,j])
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
                                               parallel=TRUE)
    k=k+1
    
  }
  
  
  return(log.scores)
  
}

