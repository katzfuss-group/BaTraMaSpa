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


# ### plotting
# pdf(file='plots/precip_samp_global.pdf',width=5.0,height=4.0)
# par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
# quilt.plot(ifelse(lon<=180,lon,lon-360),lat,x[,10],nx=288,ny=192)
# map('world',add=TRUE)
# dev.off()


## training and test data
N.test=18
N.all=ncol(x)-N.test

Ns = c(10,20,30,50,N.all)
methods=c(1,3)
splits=5
ls=array(dim=c(splits,length(Ns),length(methods)))


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
    ls[split,i.N,]=rowMeans(ls.all,na.rm=TRUE)
    
    save(ls,Ns,methods,file='output/compRes_climate_global.RData')
    
  }
  
}