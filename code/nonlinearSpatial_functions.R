
### required packages
library(invgamma)
library(fields)
library(foreach)
library(doParallel)


### parameterization functions

# nugget decay
nugfun=function(i,theta,scales) exp(theta[1]+theta[2]*log(scales[i]))

# scaling of covariates
scalingfun=function(k,theta) sqrt(exp(theta[3]*k))

# sigma (SD of nonlinear GP part)
sigmafun=function(i,theta,scales) exp(theta[4]+theta[5]*log(scales[i]))

# range
rangefun=function(theta) exp(theta[6])

# varscale (in NIG)
varscalefun=function(i,theta,scales) exp(theta[7]+theta[8]*log(scales[i]))

# concentration (in DP)
confun=function(i,theta,scales) exp(theta[9]+theta[10]*log(scales[i]))

# m as a function of theta
m.threshold=function(theta,m.max) {
  below = (scalingfun(1:m.max,theta)<.01)
  if(sum(below)==0) m=m.max else m=which.max(below)-1
  return(max(m,1))
}



## covariance kernel
kernelfun=function(X1,theta,sigma,smooth,nugget.mean=1,X2=X1){
  if(ncol(X1)==0){ return(matrix(0,nrow(X1),nrow(X2))) } else {
    # X1s = cbind(rep(mean.sd,nrow(X1)),t(t(X1)*scalingfun(1:ncol(X1),theta)))
    X1s = t(t(X1)*scalingfun(1:ncol(X1),theta))
    X2s = t(t(X2)*scalingfun(1:ncol(X2),theta))
    lin = X1s%*%t(X2s)
    nonlin= sigma^2*fields::Matern(fields::rdist(X1s,X2s),smoothness=smooth,
                                   range=rangefun(theta))
    return( (lin+nonlin)/nugget.mean )
  }
}


### compute posterior quantities
# options for mode: fit, intlik (integrated), proflik (profile)
fitMap=function(dat,NNarray.max,tuning.params,m,mode='fit',
                inds=1:nrow(dat),scales=1/(1:nrow(dat)),theta){
  
  ## tuning parameters
  if(missing(tuning.params)){
    nug.mult=4; smooth=1.5
    tuning.params=c(nug.mult,smooth)
  } else {
    nug.mult=tuning.params[1]; smooth=tuning.params[2]
  }
  
  ## extract parameters
  n=nrow(dat)
  N=ncol(dat)
  
  ## m in NNarray depends on theta
  if(missing(m)) m=m.threshold(theta,ncol(NNarray.max))
  NNarray=NNarray.max[,1:m,drop=FALSE]

  ## initialize
  if(mode=='fit'){
    chols=array(dim=c(N,N,n))
    y.tildes=array(dim=c(N,n))
    nug.means=alpha.posts=beta.posts=numeric(length=n)
  } else loglik=rep(0,n)
  
  ## loop over locations/variables
  for(i in inds){
    
    # prior of nugget
    nug.mean=nugfun(i,theta,scales)
    nugget.sd=nug.mult*nug.mean
    alpha=nug.mean^2/nugget.sd^2 + 2
    beta=nug.mean*(alpha-1)
    
    # vectors and matrices
    y = dat[i,]
    if(i==1){ K = 0 } else {
      X = t(dat[NNarray[i,1:min(i-1,m)],,drop=FALSE])
      K = kernelfun(X,theta,sigmafun(i,theta,scales),smooth,nug.mean)
    }
    G = K + diag(N)

    # expensive computations
    G.chol = tryCatch( {t(chol(G))}, error=function(cond) NULL )
    if(is.null(G.chol)){ if(mode=='fit'){stop('chol failed')} else return(-Inf) }
    y.tilde=forwardsolve(G.chol,y)
    
    # posterior of nugget
    alpha.post= alpha+N/2
    beta.post= beta + t(y.tilde)%*%y.tilde/2
    
    if(mode=='fit'){ # store quantities for later use
      chols[,,i]=G.chol; y.tildes[,i]=y.tilde; nug.means[i]=nug.mean
      alpha.posts[i]=alpha.post; beta.posts[i]=beta.post
    } else if(mode=='intlik'){ # integrated likelihood
      logdet=sum(log(diag(G.chol)))
      loglik[i]=-logdet+alpha*log(beta)-alpha.post*log(beta.post)+
        lgamma(alpha.post)-lgamma(alpha)
    } else { # profile likelihood
      nugget.hat=beta.post/(alpha.post+1) #+N/2)
      f.hat=t(y.tilde)%*%forwardsolve(G.chol,K)
      loglik[i]=sum(dnorm(y,f.hat,sqrt(nugget.hat),log=TRUE))+
        dmvnorm(f.hat,rep(0,N),K,log=TRUE)+dinvgamma(nugget,alpha,beta,log=TRUE)
    }
    
  }
  
  # return desired quantities
  if(mode=='fit'){
    return(list(chols=chols,y.tildes=y.tildes,nug.means=nug.means,
                alpha.posts=alpha.posts,beta.posts=beta.posts,scales=scales,
                dat=dat,NNarray=NNarray,theta=theta,tuning.params=tuning.params))
  } else return(sum(loglik))

}




### compute posterior quantities --- parallel version (only faster for large n/N)
# options for mode: fit, intlik (integrated), proflik (profile)
fitMapParallel=function(dat,NNarray.max,tuning.params,m,mode='fit',
                scales=1/(1:nrow(dat)),cores=NULL,theta){
  
  ## tuning parameters
  if(missing(tuning.params)){
    nug.mult=4; smooth=1.5
    tuning.params=c(nug.mult,smooth)
  } else {
    nug.mult=tuning.params[1]; smooth=tuning.params[2]
  }
  
  ## extract parameters
  n=nrow(dat)
  N=ncol(dat)
  
  ## m in NNarray depends on theta
  if(missing(m)) m=m.threshold(theta,ncol(NNarray.max))
  NNarray=NNarray.max[,1:m,drop=FALSE]
  
  ## start parallel loop over indices/locations
  registerDoParallel(cores=cores)
  par.res=foreach(i=1:n,
          .export=c('nugfun','kernelfun','sigmafun','scalingfun','rangefun')
  ) %dopar% {
    
    # prior of nugget
    nug.mean=nugfun(i,theta,scales)
    nugget.sd=nug.mult*nug.mean
    alpha=nug.mean^2/nugget.sd^2 + 2
    beta=nug.mean*(alpha-1)
    
    # vectors and matrices
    y = dat[i,]
    if(i==1){ K = 0 } else {
      X = t(dat[NNarray[i,1:min(i-1,m)],,drop=FALSE])
      K = kernelfun(X,theta,sigmafun(i,theta,scales),smooth,nug.mean)
    }
    G = K + diag(N)
    
    # expensive computations
    G.chol = tryCatch( {t(chol(G))}, error=function(cond) NULL )
    if(is.null(G.chol)){ if(mode=='fit'){stop('chol failed')} else return(-Inf) }
    y.tilde=forwardsolve(G.chol,y)
    
    # posterior of nugget
    alpha.post= alpha+N/2
    beta.post= beta + t(y.tilde)%*%y.tilde/2
    
    if(mode=='fit'){ # store quantities for later use
      list(G.chol,y.tilde,nug.mean,alpha.post,beta.post)
    } else if(mode=='intlik'){ # integrated likelihood
      logdet=sum(log(diag(G.chol)))
      loglik.i=-logdet+alpha*log(beta)-alpha.post*log(beta.post)+
        lgamma(alpha.post)-lgamma(alpha)
      loglik.i
    } else { # profile likelihood
      nugget.hat=beta.post/(alpha.post+1) #+N/2)
      f.hat=t(y.tilde)%*%forwardsolve(G.chol,K)
      loglik.i=sum(dnorm(y,f.hat,sqrt(nugget.hat),log=TRUE))+
        dmvnorm(f.hat,rep(0,N),K,log=TRUE)+dinvgamma(nugget,alpha,beta,log=TRUE)
      loglik.i
    }
    
  }
  
  
  # return desired quantities
  if(mode=='fit'){
    chols=sapply(par.res,function(x) x[[1]],simplify='array')
    y.tildes=sapply(par.res,function(x) x[[2]])
    nug.means=sapply(par.res,function(x) x[[3]])
    alpha.posts=sapply(par.res,function(x) x[[4]])
    beta.posts=sapply(par.res,function(x) x[[5]])
    return(list(chols=chols,y.tildes=y.tildes,nug.means=nug.means,
                alpha.posts=alpha.posts,beta.posts=beta.posts,scales=scales,
                dat=dat,NNarray=NNarray,theta=theta,tuning.params=tuning.params))  
  } else return(c(Reduce('+',par.res)))
  
}


### convenience function for optimization and fitting
optimFitMap=function(dat,NNarray,...,linear=FALSE,parallel=FALSE){

  ## default initial values
  theta.ini=c(log(mean(dat[1,]^2)),.2,-1,0,0,-1)
  if(linear) theta.ini=theta.ini[1:3]
  
  ## use serial or parallel function
  fitfun = if(parallel) fitMapParallel else fitMap
  
  ## integrated likelihood to optimize (on negative log scale)
  if(linear) {
    nll=function(theta) -fitfun(dat,NNarray,...,mode='intlik',
                                theta=c(theta,-Inf,0,0))
  } else {
    nll=function(theta) -fitfun(dat,NNarray,...,mode='intlik',theta=theta)
  }
  
  ## optimize
  opt=optim(theta.ini,nll,control=list(trace=0,maxit=1000,reltol=1e-3))
  
  ## theta estimate
  theta.hat = if(linear) c(opt$par,-Inf,0,0) else opt$par

  ## compute fit using theta estimate
  fit=fitfun(dat,NNarray,...,mode='fit',theta=theta.hat)
  
  return(fit)
}


### convenience function for automatic optimization and fitting
optimFitMapAuto=function(dat,locs,...){
  
  ## ordering and nearest neighbors
  ord=GPvecchia::order_maxmin_exact(locs)
  NNarray.max=GpGp::find_ordered_nn(locs[ord,],30)[,-1]
  scales=computeScales(locs[ord,],NNarray.max)
  
  ## fit posterior transport map
  fit=optimFitMap(dat[ord,],NNarray.max,scales=scales,...)
  fit$ord=ord
  
  return(fit)
}



### sampling
# options for mode: bayes or freq or score or fx or trans or invtrans or scorepm
condSamp=function(fit,mode='bayes',x.fixed=c(),obs,ind.last=nrow(fit$dat)){
  
  dat=fit$dat; NNarray=fit$NNarray; theta=fit$theta; scales=fit$scales
  nug.mult=fit$tuning.params[1]; smooth=fit$tuning.params[2]

  # extract parameters
  n=nrow(dat)
  N=ncol(dat)
  m=ncol(NNarray)
  nug.mean=fit$nug.means
  
  # loop over variables/locations
  x.new=score=c(x.fixed,rep(0,n-length(x.fixed)))
  for(i in (1+length(x.fixed)):ind.last){
    
    # predictive distribution for current sample
    if(i==1){ 
      c.star=rep(0,N) #rep(mean.sd^2,N)/nug.mean[i]
      prior.var=0
    } else {
      X = t(dat[NNarray[i,1:min(i-1,m)],,drop=FALSE])
      X.pred = 
        if(mode %in% c('score','trans','scorepm')) matrix(obs[NNarray[i,1:min(i-1,m)]],1) else 
        matrix(x.new[NNarray[i,1:min(i-1,m)]],1)
      c.star = as.numeric(kernelfun(X.pred,theta,sigmafun(i,theta,scales),
                                    smooth,nug.mean[i],X))
      prior.var=as.numeric(kernelfun(X.pred,theta,sigmafun(i,theta,scales),
                                     smooth,nug.mean[i]))
    }
    c.chol = forwardsolve(fit$chols[,,i],c.star)
    mean.pred = sum(fit$y.tildes[,i]*c.chol)
    var.pred.nonugget = prior.var - sum(c.chol^2)

    # evaluate score or sample
    if(mode=='score') {
      int.var=(fit$beta.posts[i]/fit$alpha.posts[i])*(1+var.pred.nonugget)
      score[i]=dt((obs[i]-mean.pred)/sqrt(int.var),2*fit$alpha.posts[i],log=TRUE)-.5*log(int.var)
    } else if(mode=='scorepm'){
      nugget=fit$beta.posts[i]/(fit$alpha.posts[i]-1)
      score[i]=dnorm(obs[i],mean.pred,sqrt(nugget),log=TRUE)
    } else if(mode=='fx') { 
      x.new[i]=mean.pred 
    } else if(mode=='freq'){
      nugget=fit$beta.posts[i]/(fit$alpha.posts[i]+1) #+N/2)
      x.new[i]=rnorm(1,mean.pred,sqrt(nugget))
    } else if(mode=='bayes'){
      nugget=rinvgamma(1,fit$alpha.posts[i],fit$beta.posts[i])
      x.new[i]=rnorm(1,mean.pred,sqrt(nugget*(1+var.pred.nonugget)))
    } else if(mode=='trans'){
      int.var=(fit$beta.posts[i]/fit$alpha.posts[i])*(1+var.pred.nonugget)
      x.stand=(obs[i]-mean.pred)/sqrt(int.var)
      x.new[i]=qnorm(pt(x.stand,2*fit$alpha.posts[i]))
    } else if(mode=='invtrans'){
      int.var=(fit$beta.posts[i]/fit$alpha.posts[i])*(1+var.pred.nonugget)
      x.new[i]=mean.pred+qt(pnorm(obs[i]),2*fit$alpha.posts[i])*sqrt(int.var)
    }
    
  }
  if(mode %in% c('score','scorepm')) return(sum(score[(length(x.fixed)+1):n])) else 
    return(x.new)
}


### convenience function for drawing multiple samples
NCondSamp=function(N.samp,dat,NNarray,...){
  
  fit=optimFitMap(dat,NNarray,...)
  
  x.samp=array(dim=c(nrow(dat),N.samp))
  for(j in 1:N.samp) {
    print(j)
    x.samp[,j]=condSamp(fit,mode='bayes')
  }
  
  return(x.samp)
}


### convenience function for automatic sampling
condSampAuto=function(fit,N.samp=1,...){
  
  x.samp=array(dim=c(length(fit$ord),N.samp))
  for(j in 1:N.samp) {
    x.samp[,j]=condSamp(fit,...)
  }
  
  x.samp[fit$ord,]=x.samp  # original ordering
  
  return(x.samp)
}


### plotting functions

## function for plotting and printing nonlinear parameters
plotParameters=function(theta,dat,scales=1/(1:nrow(dat))){
  par(mfrow=c(1,2))
  plot(1:n,nugfun(1:n,theta,scales),main='nug'); abline(h=var(c(dat[1,])))
  plot(1:n,sigmafun(1:n,theta,scales),main='sigma')
  par(mfrow=c(1,1))
  print(paste0('m=',m.threshold(theta,30)))
  print(paste0('range=',round(rangefun(theta),3)))
}


# ## plot nugget prior and posterior (log scale) for all indices
# plotNugget=function(fit){
# 
#   n=length(fit$nug.means)
#   yrange=log10(qinvgamma(c(.05,.95),fit$alpha.posts[c(n,1)],
#                          fit$beta.posts[c(n,1)]))
# 
#   # posterior
#   plot(log10(qinvgamma(.95,fit$alpha.posts,fit$beta.posts)),
#        xlab='index',ylab='log10(nugget)',col=3,lty=2,type='l',ylim=yrange)
#   lines(log10(qinvgamma(.5,fit$alpha.posts,fit$beta.posts)))
#   lines(log10(qinvgamma(.05,fit$alpha.posts,fit$beta.posts)),col=3,lty=2)
# 
#   # prior
#   nug.mean=fit$nug.means
#   nugget.sd=fit$tuning.params[1]*nug.mean
#   alpha=nug.mean^2/nugget.sd^2 + 2
#   beta=nug.mean*(alpha-1)
#   lines(log10(qinvgamma(.95,alpha,beta)),col=2,lty=2)
#   lines(log10(qinvgamma(.5,alpha,beta)),col=2)
#   lines(log10(qinvgamma(.05,alpha,beta)),col=2,lty=2)
# 
# }
# 
# 
# ## plot data and fitted regressions for particular i and k-th NN
# library(plotrix)
# plotGP=function(i,dat,fit,k){
#   
#   smooth=fit$tuning.params[2]
#   
#   # compute GP posterior
#   X = t(fit$dat[fit$NNarray[i,1:min(i-1,m)],,drop=FALSE])
#   K = kernelfun(X,fit$theta,sigmafun(i,fit$theta,fit$scales),
#                 smooth,nugfun(i,fit$theta,fit$scales))
#   f.hat=t(fit$y.tildes[,i])%*%forwardsolve(fit$chols[,,i],K)
#   nugget=fit$beta.posts[i]/(fit$alpha.posts[i]+1) #+N/2)
#   K.chol = forwardsolve(fit$chols[,,i],K)
#   post.var = diag(nugget*(K-t(K.chol)%*%K.chol))
#   res=fit$dat[i,]-f.hat
#   hl=1.96*sqrt(post.var)
#   
#   ## map over nearest 2 variables
#   # NN=NNarray[i,1:2]
#   # quilt.plot(fit$dat[NN[1],],fit$dat[NN[2],],fit$dat[i,])
#   
#   if(missing(k)){
#     
#     # residual plot
#     plotCI(f.hat,res,li=res-hl,ui=res+hl)
#     abline(h=0,col='grey')
#     
#   } else if(length(k)==1){
#     
#     # plot data, linear fit, and GP posterior
#     NN=fit$NNarray[i,k]
#     plot(fit$dat[NN,],fit$dat[i,],col='grey')
#     abline(lm(fit$dat[i,]~fit$dat[NN,]))
#     ord=order(fit$dat[NN,])
#     points(fit$dat[NN,ord],f.hat[ord],col=2,pch=4)
#     lines(fit$dat[NN,ord],f.hat[ord],col=2)
#     lines(fit$dat[NN,ord],f.hat[ord]-hl[ord],col=2,lty=2)
#     lines(fit$dat[NN,ord],f.hat[ord]+hl[ord],col=2,lty=2)
#     points(dat[NN,],dat[i,])
#     
#   } else {
#     
#     # make heat map over two covariates
#     NN=fit$NNarray[i,k]
#     par(mfrow=c(1,3))
#     quilt.plot(dat[NN[1],],dat[NN[2],],dat[i,])
#     quilt.plot(fit$dat[NN[1],],fit$dat[NN[2],],fit$dat[i,])
#     quilt.plot(fit$dat[NN[1],],fit$dat[NN[2],],hl)
#     par(mfrow=c(1,1))
#     
#   }
# 
# }




###########  functions for fitting model with DPM errors   ##########

#### function for cluster indicators
clust.postfun=function(clust,x,nig.par,con,output){
  
  # NIG parameters
  mu0=nig.par[1]; lambda=nig.par[2]
  alpha=nig.par[3]; beta=nig.par[4]
  
  # compute per-cluster quantities
  if(class(clust)!='factor') clust=as.factor(clust)
  K=nlevels(clust)
  nj=summary(clust)
  xj=split(x,clust)
  
  # posterior parameters for each cluster
  xbarj=sapply(xj,mean)
  xbarj[is.nan(xbarj)]=0
  ssej=sapply(xj,function(xl) sum((xl-mean(xl))^2))
  mu.post=(lambda*mu0+nj*xbarj)/(lambda+nj)
  lambda.post=lambda+nj
  alpha.post=alpha+nj/2
  beta.post=beta+.5*ssej+.5*lambda*nj*(xbarj-mu0)^2/(lambda+nj)
  
  ## output requested quantity
  if(output=='NIGsamp'){
    
    d2=invgamma::rinvgamma(K,alpha.post,beta.post)
    mus=rnorm(K,mu.post,sqrt(d2/lambda.post))
    return(rbind(mus,d2))
    
  } else {
    
    # t (marginal/predictive) params for each cluster+base
    t.df=2*c(alpha.post,alpha)
    t.mu=c(mu.post,mu0)
    t.sd=sqrt(c( beta.post/alpha.post*(lambda.post+1)/lambda.post,
                 beta/alpha*(lambda+1)/lambda ))
    
    # cluster weights
    wj=c(nj,con)/(length(clust)+con)  
    
    # output t density or clust sample
    if(output=='tfun'){
      return(function(xx) sum(wj/sum(wj)*dt((xx-t.mu)/t.sd,t.df)/t.sd))
    } else{
      t.dens=dt((output-t.mu)/t.sd,t.df)/t.sd
      levs=as.numeric(levels(clust))
      return(sample(c(levs,levs[K]+1),1,prob=wj*t.dens)) 
    }
  }
}


## update cluster indicators sequentially
update.clust=function(clust,x,nig.par,con){
  clust=as.factor(clust)
  for(i in 1:length(clust)){
    newind=clust.postfun(clust[-i],x[-i],nig.par,con,x[i])
    if(!newind %in% levels(clust)) 
      levels(clust)=c(levels(clust),newind)
    clust[i]=newind
  }
  return(as.numeric(as.factor(as.numeric(clust))))
}  

## summing on log scale
log.sum=function(x.log){
  ma=which.max(x.log)
  sum.log = x.log[ma] + log(1+sum(exp(x.log[-ma]-x.log[ma])))
  sum.log
}

## log density of normal inverse gamma
dnig=function(me,va,par) dnorm(me,par[1],sqrt(va/par[2]),log=TRUE)+
  invgamma::dinvgamma(va,par[3],par[4],log=TRUE)



### compute posterior quantities
fitDPMSerial=function(dat,NNarray.max,scales=1/(1:nrow(dat)),L=500,burnin=100,
                      thin.by=10,theta.ini,tuning.params,inds=1:nrow(dat)){
  
  ## tuning parameters
  if(missing(tuning.params)){
    nug.mult=4; smooth=1.5
    tuning.params=c(nug.mult,smooth)
  } else {
    nug.mult=tuning.params[1]; smooth=tuning.params[2]
  }
  
  ## initial values for theta parameters
  if(missing(theta.ini)) theta.ini=c(log(var(c(dat))/2),.3,-1,0,0,-.5,0,0,0,0)
  theta=theta.ini
  d=length(theta)
  
  ## extract parameters
  n=nrow(dat)
  N=ncol(dat)
  
  ## for storing thinned posterior samples
  clusts.post=array(dim=c(L/thin.by,N,n))
  clust.pars.post=array(dim=c(L/thin.by,N,2,n))
  theta.post=array(dim=c(L/thin.by,length(theta.ini)))
  
  ## initialize
  clusts=eps=array(1,dim=c(N,n))
  clust.pars=array(0,dim=c(N,2,n))
  for(i in 1:nrow(dat)) clust.pars[,2,i]=nugfun(i,theta,scales)
  gp.acc=nig.acc=con.acc=numeric(length=n)
  prop.cov=diag(.2^2,d,d)

  ## start sampler
  l.thin=1
  for(l in 1:L){
    if(l %% round(L/N*2) == 0) print(paste0('iteration ',l,' of ',L))
    
    ## propose new hyperparameter values
    theta.prop=c(theta+t(chol(prop.cov/4))%*%rnorm(length(theta)))
    m=m.threshold(theta,ncol(NNarray.max))
    m.prop=m.threshold(theta.prop,ncol(NNarray.max))
    
    ## start loop over indices/locations
    for(i in inds){
      
      y = dat[i,]
      
      ## preliminary quantities
      nug.mean=nugfun(i,theta,scales)
      nugget.sd=tuning.params[1]*nug.mean
      alpha=nug.mean^2/nugget.sd^2 + 2
      beta=nug.mean*(alpha-1)
      nig.par=c(0,varscalefun(i,theta,scales),alpha,beta)
      con=confun(i,theta,scales)
      
      ## sample eps
      if(i==1){ eps=y } else {
        mus=clust.pars[,1,i]
        D=diag(clust.pars[,2,i])
        X = t(dat[NNarray.max[i,seq_len(min(i-1,m))],,drop=FALSE])
        cov.mat=kernelfun(X,theta,sigmafun(i,theta,scales),smooth)
        G.chol=t(chol(cov.mat+D))
        GD=forwardsolve(G.chol,D)
        cov.post=D-t(GD)%*%GD
        mean.post=mus+t(GD)%*%forwardsolve(G.chol,y-mus)
        eps=c(mvtnorm::rmvnorm(1,mean.post,cov.post))
      }
      
      ## sample cluster indicators
      clusts[,i]=update.clust(clusts[,i],eps,nig.par,con)
      
      ## sample cluster parameters (NIG)
      clust.pars.unique=clust.postfun(clusts[,i],eps,nig.par,output='NIGsamp')
      clust.pars[,,i]=t(clust.pars.unique[,clusts[,i]])
      
      
      ### compute acceptance weights
      
      ## GP cov
      if(i==1){ gp.acc[i]=0 } else {
        X.prop=t(dat[NNarray.max[i,seq_len(min(i-1,m.prop))],,drop=FALSE])
        cov.mat.prop=kernelfun(X.prop,theta.prop,
                               sigmafun(i,theta.prop,scales),smooth)
        gp.acc[i]=mvtnorm::dmvnorm(y-eps,rep(0,N),cov.mat.prop,log=TRUE)-
          mvtnorm::dmvnorm(y-eps,rep(0,N),cov.mat,log=TRUE)
      }
      
      ## NIG parameters
      nug.mean=nugfun(i,theta.prop,scales)
      nugget.sd=nug.mult*nug.mean
      alpha.prop=nug.mean^2/nugget.sd^2 + 2
      beta.prop=nug.mean*(alpha.prop-1)
      nig.prop=c(0,varscalefun(i,theta.prop,scales),alpha.prop,beta.prop)
      nig.acc[i]=sum(dnig(clust.pars.unique[1,],clust.pars.unique[2,],nig.prop)-
        dnig(clust.pars.unique[1,],clust.pars.unique[2,],nig.par))
      
      ## concentration
      con.prop=confun(i,theta.prop,scales)
      K=length(unique(clusts[,i]))
      con.acc[i]=(K*log(con.prop)+lgamma(con.prop)-lgamma(N+con.prop))-
        (K*log(con)+lgamma(con)-lgamma(N+con))

    }
    
    ## accept/reject proposed hyperparameters
    gp.accept=gp.acc[!is.nan(gp.acc) & abs(gp.acc)<Inf]
    if(log(runif(1))<sum(gp.accept)) theta[3:6]=theta.prop[3:6]
    if(log(runif(1))<sum(nig.acc)) theta[c(1,2,7,8)]=theta.prop[c(1,2,7,8)]
    if(log(runif(1))<sum(con.acc)) theta[9:10]=theta.prop[9:10]

    ## save thinned iterations
    if(l%%thin.by==0){
      clusts.post[l.thin,,]=clusts
      clust.pars.post[l.thin,,,]=clust.pars
      theta.post[l.thin,]=theta
      if(l.thin>5) prop.cov=.95*2.38^2*cov(theta.post[1:l.thin,])+
                            .05*diag(.1^2,d,d)
      l.thin=l.thin+1
    }
    
  }
  
  rb=-(1:(burnin/thin.by))
  return(list(clusts.post=clusts.post[rb,,],
              clust.pars.post=clust.pars.post[rb,,,],
              theta.post=theta.post[rb,],
              dat=dat,NNarray.max=NNarray.max,
              scales=scales,tuning.params=tuning.params))
  
}




### compute posterior quantities -- IN PARALLEL
fitDPMParallel=function(dat,NNarray.max,scales=1/(1:nrow(dat)),L=500,burnin=100,
                thin.by=10,theta.ini,tuning.params,cores=NULL){
  
  ## tuning parameters
  if(missing(tuning.params)){
    nug.mult=4; smooth=1.5
    tuning.params=c(nug.mult,smooth)
  } else {
    nug.mult=tuning.params[1]; smooth=tuning.params[2]
  }
  
  ## initial values for theta parameters
  if(missing(theta.ini)) theta.ini=c(log(var(c(dat))/2),.3,-1,0,0,-.5,0,0,0,0)
  theta=theta.ini
  d=length(theta)
  
  ## extract parameters
  n=nrow(dat)
  N=ncol(dat)
  
  ## for storing thinned posterior samples
  clusts.post=array(dim=c(L/thin.by,N,n))
  clust.pars.post=array(dim=c(L/thin.by,N,2,n))
  theta.post=array(dim=c(L/thin.by,length(theta.ini)))
  
  ## initialize
  prop.cov=diag(.2^2,d,d)
  l.thin=1
  registerDoParallel(cores=cores)
  cur=list()
  for(i in 1:n) cur[[i]]=list(rep(1,N),cbind(rep(0,N),
                                             rep(nugfun(i,theta,scales),N)))
  
  ## start sampler
  for(l in 1:L){
    if(l %% round(L/N*2) == 0) print(paste0('iteration ',l,' of ',L))
    
    ## propose new hyperparameter values
    theta.prop=c(theta+t(chol(prop.cov/4))%*%rnorm(length(theta)))
    m=m.threshold(theta,ncol(NNarray.max))
    m.prop=m.threshold(theta.prop,ncol(NNarray.max))
    
    ## start loop over indices/locations
    cur=foreach(i=1:n,cur.i=iter(cur),.export=c('nugfun','varscalefun','dnig',
                  'confun','kernelfun','sigmafun','update.clust',
                  'clust.postfun','scalingfun','rangefun')
                ) %dopar% {
      
      y = dat[i,]
      
      ## preliminary quantities
      nug.mean=nugfun(i,theta,scales)
      nugget.sd=tuning.params[1]*nug.mean
      alpha=nug.mean^2/nugget.sd^2 + 2
      beta=nug.mean*(alpha-1)
      nig.par=c(0,varscalefun(i,theta,scales),alpha,beta)
      con=confun(i,theta,scales)
      
      ## sample eps
      if(i==1){ eps=y } else {
        mus=cur.i[[2]][,1]
        D=diag(cur.i[[2]][,2])
        X = t(dat[NNarray.max[i,seq_len(min(i-1,m))],,drop=FALSE])
        cov.mat=kernelfun(X,theta,sigmafun(i,theta,scales),smooth)
        G.chol=t(chol(cov.mat+D))
        G=cov.mat+D
        G.chol = tryCatch( {t(chol(G))}, error=function(cond){
          print(cond); t(chol(G+diag(1e-10*mean(diag(G)),N))) } )
        GD=forwardsolve(G.chol,D)
        cov.post=D-t(GD)%*%GD
        mean.post=mus+t(GD)%*%forwardsolve(G.chol,y-mus)
        eps=c(mvtnorm::rmvnorm(1,mean.post,cov.post))
      }
      
      ## sample cluster indicators
      clust=update.clust(cur.i[[1]],eps,nig.par,con)
      
      ## sample cluster parameters (NIG)
      clust.pars.unique=clust.postfun(clust,eps,nig.par,output='NIGsamp')
      clust.par=t(clust.pars.unique[,clust])
      
      
      ### compute acceptance weights
      
      ## GP cov
      if(i==1){ gp.acc=0 } else {
        X.prop=t(dat[NNarray.max[i,seq_len(min(i-1,m.prop))],,drop=FALSE])
        cov.mat.prop=kernelfun(X.prop,theta.prop,
                               sigmafun(i,theta.prop,scales),smooth)
        gp.acc=mvtnorm::dmvnorm(y-eps,rep(0,N),cov.mat.prop,log=TRUE)-
          mvtnorm::dmvnorm(y-eps,rep(0,N),cov.mat,log=TRUE)
        if(is.nan(gp.acc) || abs(gp.acc)==Inf) gp.acc=0
      }
      
      ## NIG parameters
      nug.mean=nugfun(i,theta.prop,scales)
      nugget.sd=nug.mult*nug.mean
      alpha.prop=nug.mean^2/nugget.sd^2 + 2
      beta.prop=nug.mean*(alpha.prop-1)
      nig.prop=c(0,varscalefun(i,theta.prop,scales),alpha.prop,beta.prop)
      nig.acc=sum(dnig(clust.pars.unique[1,],clust.pars.unique[2,],nig.prop)-
                    dnig(clust.pars.unique[1,],clust.pars.unique[2,],nig.par))
      
      ## concentration
      con.prop=confun(i,theta.prop,scales)
      K=length(unique(clust))
      con.acc=(K*log(con.prop)+lgamma(con.prop)-lgamma(N+con.prop))-
        (K*log(con)+lgamma(con)-lgamma(N+con))
      
      ## return list of results
      list(clust,clust.par,gp.acc,nig.acc,con.acc)
      
    }
    
    ## accept/reject proposed hyperparameters
    par.group=list(3:6,c(1,2,7,8),9:10)
    for(k in 1:length(par.group))
      if(log(runif(1))<sum(sapply(cur,function(x) x[[2+k]]))) 
        theta[par.group[[k]]]=theta.prop[par.group[[k]]]
    
    ## save thinned iterations
    if(l%%thin.by==0){
      clusts.post[l.thin,,]=sapply(cur,function(x) x[[1]])
      clust.pars.post[l.thin,,,]=sapply(cur,function(x) x[[2]],simplify='array')
      theta.post[l.thin,]=theta
      if(l.thin>5) prop.cov=.95*2.38^2*cov(theta.post[1:l.thin,])+
        .05*diag(.1^2,d,d)
      l.thin=l.thin+1
    }
    
  }
  
  rb=-(1:(burnin/thin.by))
  return(list(clusts.post=clusts.post[rb,,],
              clust.pars.post=clust.pars.post[rb,,,],
              theta.post=theta.post[rb,],
              dat=dat,NNarray.max=NNarray.max,
              scales=scales,tuning.params=tuning.params))
  
}


###  fitDPM: serial or parallel
fitDPM=function(...,parallel=FALSE){
  
  if(parallel) fitDPMParallel(...) else fitDPMSerial(...)
  
}




### sampling, or log score at test inputs
# options for mode: samp or score
predDPMSerial=function(fit,mode='samp',obs,x.fixed=c(),ind.last=nrow(fit$dat)){
  
  # extract quantities from fit object
  clusts.post=fit$clusts.post; clust.pars.post=fit$clust.pars.post
  theta.post=fit$theta.post; dat=fit$dat; NNarray.max=fit$NNarray.max
  scales=fit$scales; tuning.params=fit$tuning.params
  smooth=tuning.params[2]
  n=nrow(dat)
  N=ncol(dat)
  L.post=nrow(theta.post)
  
  # initialize results object
  if(mode=='samp'){
    x.samp=c(x.fixed,rep(0,n-length(x.fixed)))
  } else {
    score=rep(0,n)
    x.samp=obs
  }
  
  # loop over variables/locations
  for(i in (1+length(x.fixed)):ind.last){
    
    mix.params=array(dim=c(0,3))
    Ls=seq_len(L.post)
    if(mode=='samp') Ls=sample(Ls,1)
    for(l in Ls){
      
      # weights
      clust=factor(clusts.post[l,,i],levels=unique(clusts.post[l,,i]))
      nj=summary(clust)
      con=confun(i,theta.post[l,],scales)
      weights=c(nj,con)/(N+con)/length(Ls)
      
      # NIG params from Gibbs
      mus=clust.pars.post[l,,1,i]
      d2s=clust.pars.post[l,,2,i]
      mu=unique(mus)
      d2=unique(d2s)
      
      # draw from base
      nug.mean=nugfun(i,theta.post[l,],scales)
      nugget.sd=tuning.params[1]*nug.mean
      alpha=nug.mean^2/nugget.sd^2 + 2
      beta=nug.mean*(alpha-1)
      nig.par=c(0,varscalefun(i,theta.post[l,],scales),alpha,beta)
      d2.star=invgamma::rinvgamma(1,nig.par[3],nig.par[4])
      mu.star=rnorm(1,nig.par[1],sqrt(d2.star/nig.par[2]))
      
      # GP prediction
      m=m.threshold(theta.post[l,],ncol(NNarray.max))
      X = t(dat[NNarray.max[i,seq_len(min(i-1,m))],,drop=FALSE])
      X.pred = matrix(x.samp[NNarray.max[i,seq_len(min(i-1,m))]],1) 
      c.mat=kernelfun(X,theta.post[l,],sigmafun(i,theta.post[l,],scales),smooth)
      c.star=c(kernelfun(X.pred,theta.post[l,],sigmafun(i,theta.post[l,],scales),
                         smooth,X2=X))
      prior.var=c(kernelfun(X.pred,theta.post[l,],
                            sigmafun(i,theta.post[l,],scales),smooth))
      G.chol=t(chol( c.mat+diag(d2s) ))
      cc=forwardsolve(G.chol,c.star)
      v.star=c(prior.var-t(cc)%*%cc)
      m.star=c(t(cc)%*%forwardsolve(G.chol,dat[i,]-mus))
      
      # predictive distr: mixture mean and sd
      mm=m.star+c(mu,mu.star)
      mv=v.star+c(d2,d2.star)
      
      # mixture parameters
      mix.params=rbind(mix.params,cbind(weights,mm,mv))
      
      # # sample or log score
      # if(mode=='samp'){
      #   mi=sample(1:length(weights),1,prob=weights)
      #   x.samp[i]=rnorm(1,mm[mi],sqrt(mv[mi]))
      # } else {
      #   new.score=log.sum(log(weights)+dnorm(obs[i],mm,sqrt(mv),log=TRUE))
      #   score[i]=log.sum(c(score[i],new.score))
      # }
      
    }
    
    if(mode=='samp'){
      mi=sample(1:nrow(mix.params),1,prob=mix.params[,1])
      x.samp[i]=rnorm(1,mix.params[mi,2],sqrt(mix.params[mi,3]))
    } else
      score[i]=log.sum(log(mix.params[,1])+dnorm(obs[i],mix.params[,2],
                                                 sqrt(mix.params[,3]),log=TRUE))
    
  }
  
  if(mode=='samp') return(x.samp) else return(sum(score[(length(x.fixed)+1):n]))
  # return(mix.params)
  
}




### calculate DPM log score in parallel
predDPMParallel=function(fit,obs,cores=NULL){
  
  # extract quantities from fit object
  clusts.post=fit$clusts.post; clust.pars.post=fit$clust.pars.post
  theta.post=fit$theta.post; dat=fit$dat; NNarray.max=fit$NNarray.max
  scales=fit$scales; tuning.params=fit$tuning.params
  smooth=tuning.params[2]
  n=nrow(dat)
  N=ncol(dat)
  L.post=nrow(theta.post)
  
  # initialize results object
  score=rep(0,n)
  x.samp=obs
  
  # loop over variables/locations
  registerDoParallel(cores=cores)
  par.res=foreach(i=1:n,clusts.post.i=iapply(clusts.post,3),
                  clust.pars.post.i=iapply(clust.pars.post,4),
                  NN.i=iapply(NNarray.max,1),
                  .export=c('nugfun','varscalefun','confun','kernelfun','log.sum',
                            'sigmafun','scalingfun','rangefun','m.threshold')
  ) %dopar% {
    
    mix.params=array(dim=c(0,3))
    Ls=seq_len(L.post)
    for(l in Ls){
      
      # weights
      clust=factor(clusts.post.i[l,],levels=unique(clusts.post.i[l,]))
      nj=summary(clust)
      con=confun(i,theta.post[l,],scales)
      weights=c(nj,con)/(N+con)/length(Ls)
      
      # NIG params from Gibbs
      mus=clust.pars.post.i[l,,1]
      d2s=clust.pars.post.i[l,,2]
      mu=unique(mus)
      d2=unique(d2s)
      
      # draw from base
      nug.mean=nugfun(i,theta.post[l,],scales)
      nugget.sd=tuning.params[1]*nug.mean
      alpha=nug.mean^2/nugget.sd^2 + 2
      beta=nug.mean*(alpha-1)
      nig.par=c(0,varscalefun(i,theta.post[l,],scales),alpha,beta)
      d2.star=invgamma::rinvgamma(1,nig.par[3],nig.par[4])
      mu.star=rnorm(1,nig.par[1],sqrt(d2.star/nig.par[2]))
      
      # GP prediction
      m=m.threshold(theta.post[l,],length(NN.i))
      X = t(dat[NN.i[seq_len(min(i-1,m))],,drop=FALSE])
      X.pred = matrix(x.samp[NN.i[seq_len(min(i-1,m))]],1) 
      c.mat=kernelfun(X,theta.post[l,],sigmafun(i,theta.post[l,],scales),smooth)
      c.star=c(kernelfun(X.pred,theta.post[l,],sigmafun(i,theta.post[l,],scales),
                         smooth,X2=X))
      prior.var=c(kernelfun(X.pred,theta.post[l,],
                            sigmafun(i,theta.post[l,],scales),smooth))
      G.chol=t(chol( c.mat+diag(d2s) ))
      cc=forwardsolve(G.chol,c.star)
      v.star=c(prior.var-t(cc)%*%cc)
      m.star=c(t(cc)%*%forwardsolve(G.chol,dat[i,]-mus))
      
      # predictive distr: mixture mean and sd
      mm=m.star+c(mu,mu.star)
      mv=v.star+c(d2,d2.star)
      
      # mixture parameters
      mix.params=rbind(mix.params,cbind(weights,mm,mv))
      
    }
    
    log.sum(log(mix.params[,1])+dnorm(obs[i],mix.params[,2],
                                      sqrt(mix.params[,3]),log=TRUE))
    
  }
  
  return(c(Reduce('+',par.res)))
  
}


###  predDPM: serial or parallel
predDPM=function(...,parallel=FALSE){
  
  # note: sampling only possible in serial model (parallel=FALSE)
  if(parallel) predDPMParallel(...) else predDPMSerial(...)
  
}



##### functions for correlation-based ordering and conditioning  ####


####  compute maxmin ordering based on correlation-distance
order_maxmin_correlation <- function(cov.matrix)
{
  ## convert covariance matrix to correlation-based distance
  d=sqrt(1-abs(cov2cor(cov.matrix)))
  n=nrow(cov.matrix)
  
  ## initialize ordering
  ord=numeric(length=n)
  ord[1]=which.min(rowSums(d))
  
  ## compute maxmin ordering
  for(i in 2:n){
    if((i %% 100)==0) print(i)
    remaining=(1:n)[-ord[1:(i-1)]]
    mindist=matrixStats::rowMins(d[remaining,ord[1:(i-1)],drop=FALSE])
    # Rfast::rowMins(d[remaining,ord[1:(i-1)],drop=FALSE], value = T)
    ord[i]=remaining[which.max(mindist)]
  }
  
  return(ord)
}


#### compute m ordered NN based on correlation distance
find_nn <- function(cov.matrix,m){
  
  ## convert covariance matrix to correlation-based distance
  d=sqrt(1-abs(cov2cor(cov.matrix)))
  n=nrow(cov.matrix)
  
  ## find ordered NN
  NN=matrix(NA,n,m)
  for(i in 2:n){
    if((i %% 100)==0) print(i)
    k=min(i-1,m)
    NN[i,1:k]=order(d[i,1:(i-1)])[1:k]
  }
  
  return(NN)
}



#### compute length scales based on ordering and NNarray
computeScales=function(locs.ord,NNarray){
  n=nrow(locs.ord)
  scales=rep(0,n)
  for(i in 2:n) 
    scales[i]=rdist(locs.ord[i,,drop=FALSE],locs.ord[NNarray[i,1],,drop=FALSE])
  scales[1]=scales[2]*(scales[2]/scales[6])
  scales=scales/scales[1]
  return(scales)
}

computeScales_cor=function(cov.ord,NNarray){
  d=sqrt(1-abs(cov2cor(cov.ord)))
  n=nrow(cov.ord)
  scales=rep(0,n)
  for(i in 2:n) scales[i]=d[i,NNarray[i,1]]
  scales[1]=scales[2]*(scales[2]/scales[6])
  scales=scales/scales[1]
  return(scales)  
}


  