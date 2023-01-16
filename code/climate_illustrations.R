

###  load precip data
load(file='output/prec_anomalies.RData') # xs,locs
n=dim(xs)[1]
N=dim(xs)[2]
days=dim(xs)[3]


## correlation ordering and conditioning
m.max=30
dists=rdist(locs)
cor.est=cor(t(xs[,,1]))*exp(-dists/max(dists))
ord=order_maxmin_correlation(cor.est)
NNarray.max=find_nn(cor.est[ord,ord],m.max)
scales=computeScales_cor(cor.est[ord,ord],NNarray.max)
locs.ord=locs[ord,]
xs.ord=xs[ord,,]
dat.ord=xs.ord[,,1]



### regression scatter plot

k=2
for(i in c(40,2000)){
  pdf(file=paste0('plots/precip_','data_',i,'_3D.pdf'),width=4,height=4)
  par(mgp = c(.1,.1,0), mar=c(.1,.1,.1,.1)) # bltr
  NN=NNarray.max[i,]
  s3d=scatterplot3d::scatterplot3d(t(dat.ord[c(NN[1],NN[k],i),]),angle=30,
                                   color=rainbow(100)[cut(c(dat.ord[i,]),100)],
                                   zlab=expression(y[i]),ylab='', 
                                   pch=20,xlab='1st NN',cex.symbols = 1.2)
  s3d$box3d(col='grey')
  dims <- par("usr")
  x <- dims[1]+ 0.85*diff(dims[1:2])
  y <- dims[3]+ 0.08*diff(dims[3:4])
  text(x,y,paste0(k,'nd NN'),srt=30)
  dev.off()
}



### regression weights

weights=array(dim=c(n-2,m.max))
for(i in 3:n){
  m.i=min(i-1,m.max)
  X=t(dat.ord[NNarray.max[i,1:m.i],,drop=FALSE])
  coeff=coef(glmnet::cv.glmnet(X,dat.ord[i,],alpha=1,intercept=FALSE))
  weights[i-2,1:m.i]=coeff[-1]
}
squared.dev=apply(weights^2,2,mean,na.rm=TRUE)

pdf(file='plots/climate_sq_weights.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
plot(squared.dev[1:15],xlab='k',ylab='mean squared weights',col=1,cex=1.5)
abline(h=0)
lines(squared.dev[1:15],col=1,lwd=2)
dev.off()


### leave-one-out test data: last ensemble member
dat.train=dat.ord[,-N]
dat.test=dat.ord[,N]


### find optimal theta
fit=optimFitMap(dat.train,NNarray.max,scales=scales)
# save(fit,file='output/precip_fit.RData')
# load(file='output/precip_fit.RData')


### conditional simulation for test data
fi=c(n,300,25,0) #c(n,290,13,0)
samp.fixed=matrix(0,n,length(fi))
samp.fixed[,1]=dat.test
set.seed(1)
for(i in 2:length(fi))
  samp.fixed[,i]=condSamp(fit,x.fixed=dat.test[seq_len(fi[i])])


### plot data and conditional samples
lon=locs.ord[,1]; lat=locs.ord[,2]
maxval=max(abs(cbind(dat.train[,1:4],samp.fixed)))
ra=c(-1,1)*maxval

## plot of data
for(i in 1:4) {
  pdf(file=paste0('plots/precip_',i,'.pdf'),width=4.0,height=4.0)
  par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.1,.1)) # bltr
  quilt.plot(ifelse(lon<=180,lon,lon-360),lat,dat.train[,i],nx=37,ny=74,
             zlim=ra,add.legend=FALSE,xlab='lon',ylab='lat')
  map('world',add=TRUE)
  dev.off()
}

## plot from which to get the legend
pdf(file=paste0('plots/precip_legend.pdf'),width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.1,.1)) # bltr
quilt.plot(locs[,1],locs[,2],dat.train[,1],nx=37,ny=74,zlim=ra,
           add.legend=TRUE,xlab='lon',ylab='lat')
dev.off()

## plot conditional samples
for(i in 1:length(fi)) {
  pdf(file=paste0('plots/precip_cond_',fi[i],'.pdf'),width=4.0,height=4.0)
  par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.1,.1)) # bltr
  quilt.plot(ifelse(lon<=180,lon,lon-360),lat,samp.fixed[,i],nx=37,ny=74,
             zlim=ra,add.legend=FALSE,xlab='lon',ylab='lat')
  map('world',add=TRUE)
  dev.off()
}





##########  map/EOF coefficients  ######

### posterior median of d_i
pdf(file='plots/precip_di.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
plot(sqrt(qinvgamma(.5,fit$alpha.posts,fit$beta.posts)),
     xlab='maximin index',ylab=expression(posterior~median~of~d[i]))
dev.off()
# -> early coefficients capture more variation than late coefs


### convert to map coefficients
coef=array(dim=dim(xs))
for(year in 1:N){
  for(day in 1:days){
    print(paste0('year',year,' day',day))
    coef[,year,day]=condSamp(fit,mode='trans',obs=xs.ord[,year,day])
  }
}
# save(coef,file='output/precip_coef.RData')
# load(file='output/precip_coef.RData')

### map coefficient for held-out test field on June 1
# in original order (first by lon, then lat)
coef.test=coef[order(ord),98,1]
pdf(file='plots/precip_coef.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
plot(coef.test,xlab='index in lon/lat ordering',
     ylab='map coefficients z*')
dev.off()

acf(coef.test)


### lag-1 autocorrelation as a fct of coef index
l1ac=numeric(length=n)
for(i in 1:n) l1ac[i]=cor(c(coef[i,,2:29]),c(coef[i,,3:30]))

pdf(file='plots/precip_lag1autocor.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
plot(l1ac,xlab='coefficient index',ylab='lag-1 autocorrelation')
# abline(h=.2,v=500,col=2)
dev.off()



#### compare reconstruction to PCA

x=xs.ord[,,1]
N.test=18
N.all=ncol(x)-N.test
Ns = c(10,20,30,50,N.all)
splits=2
ls=array(dim=c(splits,length(Ns),2))
ls.k=numeric(length=N.test)

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
    
    ## fit TM
    fit=optimFitMap(data.train,NNarray.max,scales=scales)
    for(k in 1:N.test){
      x.test=data.test[,k]
      ls.k[k]=-condSamp(fit,mode='score',obs=x.test,x.fixed=x.test[1:N])
    }
    ls[split,i.N,1]=mean(ls.k)
    
    ## PCA/EOFs
    pcs=prcomp(t(data.train),tol=0)
    for(k in 1:N.test){
      x.test=data.test[,k]
      pred=predict(pcs,newdata=matrix(x.test,nrow=1))[1,]
      predfield=colSums(pred*t(pcs$rotation))
      sd.error=sqrt(mean((x.test[1:N]-predfield[1:N])^2))
      ls.k[k]=-sum(dnorm(x.test,predfield,sd.error,log=TRUE))
    }
    ls[split,i.N,2]=mean(ls.k)
    
    ## save results
    print(ls[split,i.N,])
    save(ls,Ns,file='output/reconstruction_PCA.RData')
    
  }
  
}

als=apply(ls,2:3,mean)
pdf(file='plots/recon_precip.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
matplot(Ns,als,type='l',lwd=2,xlab='n',ylab='LS',lty=1:2,col=1:2)
legend('left',c('TM','PCA'),lty=1:2,col=1:2,bg='white',lwd=2)
dev.off()





#############  UQ for subregion    ##############


## training and testing data
N.train=30
dat.train=dat.ord[,1:N.train]
dat.test=dat.ord[,(N.train+1):N]

## sample from nonlinear TM
fit=optimFitMap(dat.train,NNarray.max,scales=scales)
N.samp=500
x.samp=array(dim=c(n,N.samp))
for(j in 1:N.samp) {
  print(j)
  x.samp[,j]=condSamp(fit)
}

## sample other methods
multi.imat=autoFRK::autoFRK(Data=dat.train, loc=locs.ord, maxK= round(sqrt(n)))
covhat=multi.imat$G%*%multi.imat$M%*%t(multi.imat$G)+multi.imat$s*diag(n)
a.samp=t(mvtnorm::rmvnorm(N.samp,rep(0,n),sigma=covhat))

library(mvtnorm)
exp_nloglik <- function(params,dat){
  cov.exp=exp(params[1])*exp(-dists/exp(params[2]))
  -sum(dmvnorm(t(dat),rep(0,n),sigma=cov.exp,log=TRUE))
}
dists=rdist(locs.ord)
opt=optim(log(c(1,max(dists)*.1)),exp_nloglik,dat=dat.train,
          control=list(trace=0,maxit=100,reltol=1e-3))
covhat=exp(opt$par[1])*exp(-dists/exp(opt$par[2]))
e.samp=t(mvtnorm::rmvnorm(N.samp,rep(0,n),sigma=covhat))

## define subregion
inds=which(lon>=250 & lon<255 & lat> -13 & lat< -10)

## plot
pdf(file='plots/climate_area_fit.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
hist(colMeans(dat.test[inds,]),freq=FALSE,main='',xlab='area average',
     xlim=c(-3.5,2.5))
lines(density(colMeans(x.samp[inds,])),col=3,lwd=2,lty=1)
lines(density(colMeans(a.samp[inds,])),col=2,lwd=2,lty=2)
lines(density(colMeans(e.samp[inds,])),col=4,lwd=2,lty=3)
legend('topleft',c('nonlin','expCov','autoFRK'),col=c(3,2,4),lty=c(1,2,3))
dev.off()
