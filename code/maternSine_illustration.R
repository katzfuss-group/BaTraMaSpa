
#####  illustration of the simulated matern+sine model  #####


### create the matern+sine model and sample from it
n=30^2
Ns=100
N.test=1000
nonlin=2 # degree of nonlinearity
reg=TRUE # regular grid or randomly sampled locations
sepa=0 # separation for bimodal residual errors
set.seed(99)
source('code/createMaternSine.R')



#### plot one realization of field
pdf(file=paste0('plots/sine_','field','.pdf'),width=5,height=4.0)
par(mgp = c(1.6,.5,.5), mar=c(.1,.1,.1,.3)) # bltr
quilt.plot(locs.ord[,1],locs.ord[,2],data.all[,1],
           nx=sqrt(n),ny=sqrt(n),xaxt='n',yaxt='n')
dev.off()




#### train linear and nonlinear models

N=100
dat.train=data.all[,1:N]
scales=computeScales(locs.ord,NNarray.max)

fit.lin=optimFitMap(dat.train,NNarray.max,scales=scales,linear=TRUE)
fit.nonlin=optimFitMap(dat.train,NNarray.max,scales=scales,linear=FALSE)
# save(fit.lin,fit.nonlin,file='output/sine_singlefit.RData')
# load(file='output/sine_singlefit.RData')
  



#### plot fit for a single location

i=80

## obtain prediction/fct surface
NN=as.numeric(na.omit(NNarray.max[i,]))
x.fixed=rep(0,i-1)
x.fixed[NN]=rowMeans(dat.train[NN,])
n.vals=20
NN.vals=array(dim=c(2,n.vals))
for(j in 1:2) {
  qs=range(dat.train[NN[j],]) #quantile(data.all[NN[j],],c(.05,.95))
  NN.vals[j,]=seq(qs[1],qs[2],length=n.vals)
}
fx=fx.lin=true.fun=array(dim=c(n.vals,n.vals))
for(k in 1:n.vals){
  for(l in 1:n.vals){
    x.fixed[NN[1:2]]=c(NN.vals[1,k],NN.vals[2,l])
    fx[k,l]=condSamp(fit.nonlin,mode='fx',x.fixed,ind.last=i)[i]
    fx.lin[k,l]=condSamp(fit.lin,mode='fx',x.fixed,ind.last=i)[i]
    true.fun[k,l]=f(x.fixed[NN],weights[i])
  }
}

  
## remove unobserved areas
temp=quilt.plot(dat.train[NN[1],],dat.train[NN[2],],dat.train[i,],
                nx=n.vals,ny=n.vals,plot=FALSE)$z
wi=numeric(length=n.vals)
for(l in 1:n.vals) wi[l]=max(abs(range(which(!is.na(temp[,l])))-l))
bwi=max(wi[abs(wi)!=Inf])
for(k in 1:n.vals){
  for(l in 1:n.vals){
    if(abs(k-l)>bwi) fx[k,l]=fx.lin[k,l]=true.fun[k,l]=NA
  }
}


## save plots

ra=range(c(range(cbind(true.fun,fx,fx.lin),na.rm=TRUE),range(dat.train[i,])))

pdf(file=paste0('plots/sine_bw_','data','.pdf'),width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.1,.1)) # bltr
quilt.plot(dat.train[NN[1],],dat.train[NN[2],],dat.train[i,],
           zlim=ra,add.legend=FALSE,
           xlab='1st NN',ylab='2nd NN',nx=n.vals,ny=n.vals)
dev.off()

plotnames=c('true','nonlin','lin')
plotvals=list(true.fun,fx,fx.lin)
for(i in 1:3) {
  pdf(file=paste0('plots/sine_bw_',plotnames[i],'.pdf'),width=4.0,height=4.0)
  par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.1,.1)) # bltr
  image(NN.vals[1,],NN.vals[2,],plotvals[[i]],zlim=ra,col=tim.colors(64),
        xlab='1st NN',ylab='2nd NN')
  dev.off()
}

# plot from which to get the legend
pdf(file=paste0('plots/sine_bw_legend.pdf'),width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(0,0,0,0)) # bltr
image.plot(NN.vals[1,],NN.vals[2,],plotvals[[i]],zlim=ra,legend.only=TRUE)
dev.off()


# ### version using 3d scatterplot
# library(scatterplot3d)
# k=2
# s3d=scatterplot3d(t(dat.train[c(NN[1],NN[k],i),]),angle=30,
#                   xlim=ra,ylim=ra,zlim=ra,pch=20,
#                   color=rainbow(100)[cut(dat.train[i,],100)],
#                   xlab='1st NN',ylab=paste0(k,'th NN'),zlab=expression(y[i]))
# s3d$box3d(col='grey')
# 
# for(l in 1:3) {
#   s3d=scatterplot3d(cbind(expand.grid(NN.vals[1,],NN.vals[2,]),c(plotvals[[l]])),
#                     angle=30,xlim=ra,ylim=ra,zlim=ra,
#                     color=rainbow(100)[cut(c(plotvals[[l]]),100)],
#                     xlab='1st NN',ylab=paste0(k,'th NN'),zlab=expression(y[i]))
#   s3d$box3d(col='grey')
# }


## m for nonlinear fit
m.threshold(fit.nonlin$theta,m.max)



###### transformation to latent space

## exact map
z=rep(0,n)
for(i in 1:n){
  fx = f(y[NNarray.max[i,]],weights[i,])
  z[i] = (test.data[i,j]-fx)/cond.sds[i]
}

## transformation
data.latent=test.data=data.test
for(j in 1:ncol(test.data)) { print(j)
  data.latent[,j]=condSamp(fit.nonlin,mode='trans',obs=test.data[,j])
}
inds=seq(1,N.test,by=2)
avg.orig=avg.latent=array(dim=c(n,length(inds)))
for(k in 1:length(inds)){
  indsk=if(k<length(inds)) inds[k]:(inds[k+1]-1) else inds[k]:N.test
  avg.orig[,k]=rowMeans(test.data[,indsk])
  avg.latent[,k]=condSamp(fit.nonlin,mode='invtrans',
                          obs=rowMeans(data.latent[,indsk]))
}
# save(data.latent,avg.orig,avg.latent,file='output/materndat_trans.RData')
# load(file='output/materndat_trans.RData')


### plot transformed data
pdf(file='plots/matern_trans_data.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
plot(data.latent[,10],xlab='index i in maximin ordering',
     ylab='transformed data')
dev.off()

pdf(file='plots/matern_trans_qq.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
qqnorm(data.latent[,10],main='',xlab='Theoretical Gaussian Quantiles')
dev.off()



### vs NN

i=80
NN=as.numeric(NNarray.max[i,1:2])
ras=apply(rbind(apply(test.data[c(NN,i),],1,range),
                apply(avg.orig[c(NN,i),],1,range),
                apply(avg.latent[c(NN,i),],1,range)),2,range)

pdf(file=paste0('plots/sine_avg_test.pdf'),width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.1,.1)) # bltr
quilt.plot(test.data[NN[1],],test.data[NN[2],],test.data[i,],add.legend=FALSE,
           xlim=ras[,1],ylim=ras[,2],zlim=ras[,3],
           xlab='1st NN',ylab='2nd NN')
dev.off()

pdf(file=paste0('plots/sine_avg_orig.pdf'),width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.1,.1)) # bltr
quilt.plot(avg.orig[NN[1],],avg.orig[NN[2],],avg.orig[i,],add.legend=FALSE,
           xlim=ras[,1],ylim=ras[,2],zlim=ras[,3],
           xlab='1st NN',ylab='2nd NN')
dev.off()

pdf(file=paste0('plots/sine_avg_latent.pdf'),width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.1,.1)) # bltr
quilt.plot(avg.latent[NN[1],],avg.latent[NN[2],],avg.latent[i,],add.legend=FALSE,
           xlim=ras[,1],ylim=ras[,2],zlim=ras[,3],
           xlab='1st NN',ylab='2nd NN')
dev.off()

# plot from which to get the legend
pdf(file=paste0('plots/sine_avg_legend.pdf'),width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(0,0,0,0)) # bltr
image.plot(avg.latent[NN[1],],avg.latent[NN[2],],avg.latent[i,],
           zlim=ras[,3],legend.only=TRUE)
dev.off()





### plot realization of matern+sine+bimodal
n=30^2
Ns=2
N.test=2
nonlin=1 # linear or nonlinear?
reg=TRUE # regular grid or randomly sampled locations
sepa=3.5 # separation for bimodal residual errors
set.seed(99)
source('code/createMaternSine.R')

pdf(file=paste0('plots/sine_','field_bimodal','.pdf'),width=5,height=4.0)
par(mgp = c(1.6,.5,.5), mar=c(.1,.1,.1,.3)) # bltr
quilt.plot(locs.ord[,1],locs.ord[,2],data.all[,1],
           nx=sqrt(n),ny=sqrt(n),xaxt='n',yaxt='n')
dev.off()
