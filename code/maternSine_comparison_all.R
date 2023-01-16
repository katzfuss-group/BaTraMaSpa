#####  nonlinear sine map based on matern   #####

# setwd("G:/My Drive/projects/triangular_spatial")
source('code/nonlinearSpatial_functions.R')
source('code/logScore_comparison.R')



############   small, regular setting: linear

###  create maternSine model
n=30^2
Ns=c(5,10,20,30,50,100,200)
N.test=50
nonlin=0 # degree of nonlinearity
reg=TRUE # regular grid or randomly sampled locations
sepa=0 # separation for bimodal residual errors
source('code/createMaternSine.R')


### plot one realization of field
pdf(file='plots/simLR900.pdf',width=5,height=4.0)
par(mgp = c(1.6,.5,.5), mar=c(.1,.1,.1,.3)) # bltr
quilt.plot(locs.ord[,1],locs.ord[,2],data.all[,1],
           nx=sqrt(n),ny=sqrt(n),xaxt='n',yaxt='n')
dev.off()

## taper matrix
dists=rdist(locs.ord)
tap=exp(-dists/max(dists))

## grid and windows for local matern method
grdwin=local_windows(locs)
ro=order(ord)

### compute log scores for different methods
ls.exact=numeric(length=N.test)
methods=1:4
kls=array(dim=c(length(Ns),length(methods)+4))
for(i.N in 1:length(Ns)){
  
  ## training data
  N=Ns[i.N]; print(paste0('N=',N))
  data.train=data.all[,1:N]
  
  ## LS for transport map methods
  ls=compScore(locs.ord,data.train,data.test,methods=methods,
               NNarray.max=NNarray.max)
  
  ## LS for matern GP
  ls=rbind(ls,matern_scores(locs.ord,data.train,data.test))
  
  ## local matern
  ls=rbind(ls,local_scores(grdwin,data.train[ro,],data.test[ro,]))
  
  ## tapered sample covariance
  cov.est=cov(t(data.train))*tap
  ls=rbind(ls,mean(dmvnorm(t(data.test),rep(0,n),sigma=cov.est,log=TRUE)))

  ## autoFRK
  multi.imat=autoFRK::autoFRK(Data=data.train,loc=locs.ord,maxK= round(sqrt(n)))
  covhat=multi.imat$G%*%multi.imat$M%*%t(multi.imat$G)+multi.imat$s*diag(n)
  ls=rbind(ls,mean(dmvnorm(t(data.test),rep(0,n),sigma=covhat,log=TRUE)))
  
  ## compute compute KLs as difference in avg log scores
  for(j in 1:N.test) ls.exact[j]=score.exact(data.test[,j])
  kls[i.N,]=mean(ls.exact)-rowMeans(ls,na.rm=TRUE)
  
  save(kls,Ns,methods,file='output/compRes_all_sin_lin.RData')
  
}



############   small, regular setting: nonlinear

###  create maternSine model
n=30^2
Ns=c(5,10,20,30,50,100,200)
N.test=50
nonlin=2 # degree of nonlinearity
reg=TRUE # regular grid or randomly sampled locations
sepa=0 # separation for bimodal residual errors
source('code/createMaternSine.R')


### plot one realization of field
pdf(file='plots/simNR900.pdf',width=5,height=4.0)
par(mgp = c(1.6,.5,.5), mar=c(.1,.1,.1,.3)) # bltr
quilt.plot(locs.ord[,1],locs.ord[,2],data.all[,1],
           nx=sqrt(n),ny=sqrt(n),xaxt='n',yaxt='n')
dev.off()


## taper matrix
dists=rdist(locs.ord)
tap=exp(-dists/max(dists))

## grid and windows for local matern method
grdwin=local_windows(locs)
ro=order(ord)

### compute log scores for different methods
ls.exact=numeric(length=N.test)
methods=1:4
kls=array(dim=c(length(Ns),length(methods)+4))
for(i.N in 1:length(Ns)){
  
  ## training data
  N=Ns[i.N]; print(paste0('N=',N))
  data.train=data.all[,1:N]
  
  ## LS for transport map methods
  ls=compScore(locs.ord,data.train,data.test,methods=methods,
               NNarray.max=NNarray.max)
  
  ## LS for matern GP
  ls=rbind(ls,matern_scores(locs.ord,data.train,data.test))
  
  ## local matern
  ls=rbind(ls,local_scores(grdwin,data.train[ro,],data.test[ro,]))
  
  ## tapered sample covariance
  cov.est=cov(t(data.train))*tap
  ls=rbind(ls,mean(dmvnorm(t(data.test),rep(0,n),sigma=cov.est,log=TRUE)))
  
  ## autoFRK
  multi.imat=autoFRK::autoFRK(Data=data.train,loc=locs.ord,maxK= round(sqrt(n)))
  covhat=multi.imat$G%*%multi.imat$M%*%t(multi.imat$G)+multi.imat$s*diag(n)
  ls=rbind(ls,mean(dmvnorm(t(data.test),rep(0,n),sigma=covhat,log=TRUE)))
  
  ## compute compute KLs as difference in avg log scores
  for(j in 1:N.test) ls.exact[j]=score.exact(data.test[,j])
  kls[i.N,]=mean(ls.exact)-rowMeans(ls,na.rm=TRUE)
  
  save(kls,Ns,methods,file='output/compRes_all_sin.RData')
  
}




################  large, irregular setting

###  create maternSine model
n=80^2
Ns=c(5,10,20,30,50,100,200)
N.test=50
nonlin=2 # degree of nonlinearity
reg=FALSE # regular grid or randomly sampled locations
sepa=0 # separation for bimodal residual errors
source('code/createMaternSine.R')


### plot one realization of field
pdf(file='plots/simNI3600.pdf',width=5,height=4.0)
par(mgp = c(1.6,.5,.5), mar=c(.1,.1,.1,.3)) # bltr
quilt.plot(locs.ord[,1],locs.ord[,2],data.all[,1],
           nx=200,ny=200,xaxt='n',yaxt='n')
dev.off()

## taper matrix
dists=rdist(locs.ord)
tap=exp(-dists/max(dists))


### compute log scores for different methods
ls.exact=numeric(length=N.test)
methods=1:4
kls=array(dim=c(length(Ns),length(methods)+3))
for(i.N in 1:length(Ns)){
  
  ## training data
  N=Ns[i.N]; print(paste0('N=',N))
  data.train=data.all[,1:N]
  
  ## LS for transport map methods
  ls=compScore(locs.ord,data.train,data.test,methods=methods,NNarray.max)
  
  ## LS for matern GP
  ls=rbind(ls,matern_scores(locs.ord,data.train,data.test))
  
  ## tapered sample covariance
  cov.est=cov(t(data.train))*tap
  ls=rbind(ls,mean(dmvnorm(t(data.test),rep(0,n),sigma=cov.est,log=TRUE)))
  
  ## autoFRK
  multi.imat=autoFRK::autoFRK(Data=data.train,loc=locs.ord,maxK= round(sqrt(n)))
  covhat=multi.imat$G%*%multi.imat$M%*%t(multi.imat$G)+multi.imat$s*diag(n)
  ls=rbind(ls,mean(dmvnorm(t(data.test),rep(0,n),sigma=covhat,log=TRUE)))
  
  ## compute compute KLs as difference in avg log scores
  for(j in 1:N.test) ls.exact[j]=score.exact(data.test[,j])
  kls[i.N,]=mean(ls.exact)-rowMeans(ls,na.rm=TRUE)
  
  save(kls,Ns,methods,n,reg,file='output/compRes_all_sin_random.RData')
  
}



############   small, regular, **bimodal**

###  create maternSine model
n=30^2
Ns=c(5,10,20,30,50,100,200)
N.test=50
nonlin=1 # degree of nonlinearity
reg=TRUE # regular grid or randomly sampled locations
sepa=3.5 # separation for bimodal residual errors
source('code/createMaternSine.R')


### plot one realization of field
pdf(file='plots/simNR900B.pdf',width=5,height=4.0)
par(mgp = c(1.6,.5,.5), mar=c(.1,.1,.1,.3)) # bltr
quilt.plot(locs.ord[,1],locs.ord[,2],data.all[,2],
           nx=sqrt(n),ny=sqrt(n),xaxt='n',yaxt='n')
dev.off()

## taper matrix
dists=rdist(locs.ord)
tap=exp(-dists/max(dists))

## grid and windows for local matern method
grdwin=local_windows(locs)
ro=order(ord)

### compute log scores for different methods
ls.exact=numeric(length=N.test)
methods=1:5
kls=array(dim=c(length(Ns),length(methods)+4))
for(i.N in 1:length(Ns)){
  
  ## training data
  N=Ns[i.N]; print(paste0('N=',N))
  data.train=data.all[,1:N]
  
  ## LS for transport map methods
  ls=compScore(locs.ord,data.train,data.test,methods=methods,NNarray.max)
  
  ## LS for matern GP
  ls=rbind(ls,matern_scores(locs.ord,data.train,data.test))
  
  ## local matern
  ls=rbind(ls,local_scores(grdwin,data.train[ro,],data.test[ro,]))
  
  ## tapered sample covariance
  cov.est=cov(t(data.train))*tap
  ls=rbind(ls,mean(dmvnorm(t(data.test),rep(0,n),sigma=cov.est,log=TRUE)))
  
  ## autoFRK
  multi.imat=autoFRK::autoFRK(Data=data.train,loc=locs.ord,maxK= round(sqrt(n)))
  covhat=multi.imat$G%*%multi.imat$M%*%t(multi.imat$G)+multi.imat$s*diag(n)
  ls=rbind(ls,mean(dmvnorm(t(data.test),rep(0,n),sigma=covhat,log=TRUE)))
  
  ## compute compute KLs as difference in avg log scores
  for(j in 1:N.test) ls.exact[j]=score.exact(data.test[,j])
  kls[i.N,]=mean(ls.exact)-rowMeans(ls,na.rm=TRUE)
  
  save(kls,Ns,methods,file='output/compRes_all_sin_DPM.RData')
  
}



############   another non-Gaussian setting

n=30^2
Ns=c(5,10,20,30,50,100,200)
N.test=50
N.total=max(Ns)+N.test

grid.oneside=seq(0,1,length.out=sqrt(n))
locs=as.matrix(expand.grid(grid.oneside,grid.oneside))
dist=fields::rdist(locs)

range=.3; smooth=.5
cov1=fields::Matern(dist,range=range,nu=smooth)

y1=t(chol(cov1))%*%matrix(rnorm(n*N.total),nrow=n) 
y2=matrix(sin(apply(locs-.5,1,sum)/.05),nrow=n,ncol=N.total)
data.all=y1*y2
data.test=data.all[,-(1:max(Ns))]


### plot one realization of field
pdf(file='plots/sim_prod.pdf',width=5,height=4.0)
par(mgp = c(1.6,.5,.5), mar=c(.1,.1,.1,.3)) # bltr
quilt.plot(locs[,1],locs[,2],data.all[,1],
           nx=sqrt(n),ny=sqrt(n),xaxt='n',yaxt='n')
dev.off()

## taper matrix
tap=exp(-dist/max(dist))

## grid and windows for local matern method
grdwin=local_windows(locs)

### compute log scores for different methods
methods=1:4
ls=array(dim=c(length(Ns),length(methods)+4))
for(i.N in 1:length(Ns)){
  
  ## training data
  N=Ns[i.N]; print(paste0('N=',N))
  data.train=data.all[,1:N]
  
  ## compute avg log scores
  ls.all=compScore(locs,data.train,data.test,methods=methods)
  ls[i.N,methods]=-rowMeans(ls.all,na.rm=TRUE)
  
  ## matern covariance
  ls[i.N,length(methods)+1]=-matern_scores(locs,data.train,data.test)
  
  ## local matern
  ls[i.N,length(methods)+2]=-mean(local_scores(grdwin,data.train,data.test))
  
  ## tapered sample covariance
  cov.est=cov(t(data.train))*tap
  ls[i.N,length(methods)+3]=-mean(dmvnorm(t(data.test),rep(0,n),sigma=cov.est,
                                          log=TRUE))
  
  ## autoFRK
  multi.imat=autoFRK::autoFRK(Data=data.train,loc=locs,maxK= round(sqrt(n)))
  covhat=multi.imat$G%*%multi.imat$M%*%t(multi.imat$G)+multi.imat$s*diag(n)
  ls[i.N,length(methods)+4]=-mean(dmvnorm(t(data.test),rep(0,n),sigma=covhat,
                                          log=TRUE))
  
  save(ls,Ns,methods,file='output/compRes_all_prod.RData')
  
}
