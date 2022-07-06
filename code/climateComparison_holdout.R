#####  comparison on climate data   #####

# setwd("G:/My Drive/projects/triangular_spatial")
source('code/nonlinearSpatial_functions.R')
source('code/logScore_comparison.R')



## load climate data: precipitation anomalies july 1
load(file='output/prec_anomalies.RData') # xs,locs
x=xs[,,1] # 1st day of july
n=nrow(x)

## grid and windows for local matern method
grdwin=local_windows(locs)

## 3D coordinates for chordal (Euclidean) distance
l=locs[,1]/360*2*pi; L=locs[,2]/360*2*pi
locs=cbind(cos(L)*cos(l),cos(L)*sin(l),sin(L))

## taper matrix
dists=rdist(locs)
tap=exp(-dists/max(dists))


## training and test data
N.test=18
N.all=ncol(x)-N.test

Ns = c(10,20,30,50,N.all)
methods=1:5
splits=5
ls=array(dim=c(splits,length(Ns),9,2))
cor.ord=FALSE

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
    s.all=compScore(locs,data.train,data.test,methods=methods,cor.ord=cor.ord,
                    holdout=TRUE)
    ls[split,i.N,methods,1]=rowMeans(s.all$log.scores,na.rm=TRUE)
    ls[split,i.N,methods,2]=rowMeans(s.all$ho.scores,na.rm=TRUE)
    
    ## matern covariance
    ls[split,i.N,6,]=matern_scores(locs,data.train,data.test,holdout=s.all$hoi)
    
    ## local matern
    ls[split,i.N,7,]=local_scores(grdwin,data.train,data.test,holdout=s.all$hoi)
    
    ## tapered sample covariance
    cov.est=cov(t(data.train))*tap
    ls[split,i.N,8,1]=mean(dmvnorm(t(data.test),rep(0,n),sigma=cov.est,log=TRUE))
    ls[split,i.N,8,2]=ls.ho.gaussian(cov.est,data.test,holdout=s.all$hoi)
    
    ## autoFRK
    multi.imat=autoFRK::autoFRK(Data=data.train, loc=locs, maxK= round(sqrt(n)))
    covhat=multi.imat$G%*%multi.imat$M%*%t(multi.imat$G)+multi.imat$s*diag(n)
    ls[split,i.N,9,1]=mean(dmvnorm(t(data.test),rep(0,n),sigma=covhat,log=TRUE))
    ls[split,i.N,9,2]=ls.ho.gaussian(covhat,data.test,holdout=s.all$hoi)
    
    ## save results
    print(ls[split,i.N,,])
    save(ls,Ns,file='output/compRes_climate_holdout.RData')
    
  }
  
}


