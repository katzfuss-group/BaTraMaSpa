#####  comparison on climate data   #####

# setwd("G:/My Drive/projects/triangular_spatial")
source('code/nonlinearSpatial_functions.R')
source('code/logScore_comparison.R')



## load climate data: precipitation anomalies july 1
load(file='output/prec_anomalies.RData') # xs,locs
x=xs[,,1] # 1st day of july
n=nrow(x)

## taper matrix
dists=rdist(locs)
tap=exp(-dists/max(dists))

## likelihood for GP with exponential cov
library(mvtnorm)
exp_nloglik <- function(params,dat){
  cov.exp=exp(params[1])*exp(-dists/exp(params[2]))
  -sum(dmvnorm(t(dat),rep(0,n),sigma=cov.exp,log=TRUE))
}

## training and test data
N.test=18
N.all=ncol(x)-N.test

Ns = c(10,20,30,50,N.all)
methods=1:5
splits=5
ls=array(dim=c(splits,length(Ns),8))
cor.ord=TRUE

## random splits
for(split in 1:splits){
  
  perm=sample(1:ncol(x),ncol(x))
  data.all=x[,perm[1:N.all]]
  data.test=x[,perm[(N.all+1):length(perm)]]
  
  ## compute log scores for different methods
  for(i.N in 1:length(Ns)){
    
    ## training data
    N=Ns[i.N]; print(paste0('N=',N))
    data.train=data.all[,1:N]
    
    ## compute avg log scores
    ls.all=compScore(locs,data.train,data.test,methods=methods,cor.ord=cor.ord)
    ls[split,i.N,methods]=rowMeans(ls.all,na.rm=TRUE)
    
    ## tapered sample covariance
    cov.est=cov(t(data.train))*tap
    ls[split,i.N,6]=sum(dmvnorm(t(data.test),rep(0,n),
                                 sigma=cov.est,log=TRUE))
    
    ## exponential covariance
    opt=optim(log(c(1,max(dists)*.1)),exp_nloglik,dat=data.train,
              control=list(trace=0,maxit=100,reltol=1e-3))
    par.est=opt$par
    ls[split,i.N,7]=-exp_nloglik(par.est,data.test)
    
    ## autoFRK
    multi.imat=autoFRK::autoFRK(Data=data.train, loc=locs, maxK= round(sqrt(n)))
    covhat=multi.imat$G%*%multi.imat$M%*%t(multi.imat$G)+multi.imat$s*diag(n)
    ls[split,i.N,8]=sum(dmvnorm(t(data.test),rep(0,n),sigma=covhat,log=TRUE))
    
    ## save results
    save(ls,Ns,file='output/compRes_climate_scales.RData')
    
  }
  
}


