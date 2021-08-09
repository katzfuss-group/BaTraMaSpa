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


### compute log scores for different methods
ls.exact=numeric(length=N.test)
methods=1:4
kls=array(dim=c(length(Ns),length(methods)))
for(i.N in 1:length(Ns)){
  
  ## training data
  N=Ns[i.N]; print(paste0('N=',N))
  data.train=data.all[,1:N]
  
  ## compute compute KLs as difference in avg log scores
  ls=compScore(locs.ord,data.train,data.test,methods=methods,NNarray.max)
  for(j in 1:N.test) ls.exact[j]=score.exact(data.test[,j])
  kls[i.N,]=mean(ls.exact)-rowMeans(ls,na.rm=TRUE)
  
  save(kls,Ns,methods,file='output/compRes_sin_lin.RData')
  
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


### compute log scores for different methods
ls.exact=numeric(length=N.test)
methods=1:4
kls=array(dim=c(length(Ns),length(methods)))
for(i.N in 1:length(Ns)){
  
  ## training data
  N=Ns[i.N]; print(paste0('N=',N))
  data.train=data.all[,1:N]
  
  ## compute compute KLs as difference in avg log scores
  ls=compScore(locs.ord,data.train,data.test,methods=methods,NNarray.max)
  for(j in 1:N.test) ls.exact[j]=score.exact(data.test[,j])
  kls[i.N,]=mean(ls.exact)-rowMeans(ls,na.rm=TRUE)
  
  save(kls,Ns,methods,file='output/compRes_sin.RData')
  
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


### compute log scores for different methods
ls.exact=numeric(length=N.test)
methods=1:4
kls=array(dim=c(length(Ns),length(methods)))
for(i.N in 1:length(Ns)){
  
  ## training data
  N=Ns[i.N]; print(paste0('N=',N))
  data.train=data.all[,1:N]
  
  ## compute compute KLs as difference in avg log scores
  ls=compScore(locs.ord,data.train,data.test,methods=methods,NNarray.max)
  for(j in 1:N.test) ls.exact[j]=score.exact(data.test[,j])
  kls[i.N,]=mean(ls.exact)-rowMeans(ls,na.rm=TRUE)
  
  save(kls,Ns,methods,n,reg,file='output/compRes_sin_random.RData')
  
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


### compute log scores for different methods
ls.exact=numeric(length=N.test)
methods=1:5
kls=array(dim=c(length(Ns),length(methods)))
for(i.N in 1:length(Ns)){
  
  ## training data
  N=Ns[i.N]; print(paste0('N=',N))
  data.train=data.all[,1:N]
  
  ## compute compute KLs as difference in avg log scores
  ls=compScore(locs.ord,data.train,data.test,methods=methods,NNarray.max)
  for(j in 1:N.test) ls.exact[j]=score.exact(data.test[,j])
  kls[i.N,]=mean(ls.exact)-rowMeans(ls,na.rm=TRUE)
  
  save(kls,Ns,methods,file='output/compRes_sin_DPM.RData')
  
}

