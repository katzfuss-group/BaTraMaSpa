################   global analysis   ##################

# setwd("G:/My Drive/projects/triangular_spatial")
source('code/nonlinearSpatial_functions.R')
source('code/logScore_comparison.R')


#### load climate data: precipitation rate july 1

## global data
load(file='data/prec_all.RData') # lon,lat,prec
x=log(prec+1e-10)
x=t(scale(t(x))) # standardize

## order & condition based on 3D coordinates
l=lon/360*2*pi; L=lat/360*2*pi
locs=cbind(cos(L)*cos(l),cos(L)*sin(l),sin(L))
cor.ord=FALSE


## training and test data
N.test=18
N.all=ncol(x)-N.test

Ns = c(10,20,30,50,N.all)
methods=c(1,3)
splits=5

split=1

  perm=sample(1:ncol(x),ncol(x))
  data.all=x[,perm[1:N.all]]
  data.test=x[,perm[(N.all+1):length(perm)]]
  
  i.N=length(Ns)

    ## training data
    N=Ns[i.N]; print(paste0('N=',N))
    data.train=data.all[,1:N]
    
    
m.max=30
      
      ## general settings and parameters
      n=nrow(data.train)
      N=ncol(data.train)
      N.test=ncol(data.test)
      
          ord=GPvecchia::order_maxmin_exact(locs)
          NNarray.max=GpGp::find_ordered_nn(locs[ord,],m.max)[,-1]
          scales=computeScales(locs[ord,],NNarray.max)
        dat.train=data.train[ord,]
        dat.test=data.test[ord,]

      
        # train model
        fit=optimFitMap(dat.train,NNarray.max,scales=scales,
                        linear=FALSE,parallel=FALSE)

fg=list(theta=fit$theta,alpha.posts=fit$alpha.posts,beta.posts=fit$beta.posts)        
save(n,N,fg,file='output/fit_global.RData')
