

#####  plot global precip data   #####

# setwd("G:/My Drive/projects/triangular_spatial")
source('code/nonlinearSpatial_functions.R')


#### load climate data: precipitation rate july 1

## global data
load(file='data/prec_all.RData') # lon,lat,prec
x=log(prec+1e-10)
x=t(scale(t(x))) # standardize

long=ifelse(lon<=180,lon,lon-360)


####### conditional simulation for test data #####


### fit map to first 97 samples
l=lon/360*2*pi; L=lat/360*2*pi
locs=cbind(cos(L)*cos(l),cos(L)*sin(l),sin(L))

ord=GPvecchia::order_maxmin_exact(locs)
NNarray.max=GpGp::find_ordered_nn(locs[ord,],30)[,-1]
scales=computeScales(locs[ord,],NNarray.max)

dat.train=x[ord,-ncol(x)]
dat.test=x[ord,ncol(x)]

fit=optimFitMap(dat.train,NNarray.max,scales=scales)
# save(long,lat,ord,dat.test,fit,file='output/climate_global_fit97.RData')
# load(file='output/climate_global_fit97.RData')


### conditional simulation for 98th sample
n=length(dat.test)
fi=c(n,5000,500,0)
samp.fixed=matrix(0,n,length(fi))
samp.fixed[,1]=dat.test
set.seed(1)
for(i in 2:length(fi))
  samp.fixed[,i]=condSamp(fit,x.fixed=dat.test[seq_len(fi[i])])


### plot data and conditional samples
library(oce)
data(coastlineWorld)
maxval=quantile(abs(cbind(x[,1:4],samp.fixed)),.999)
cm=colormap(breaks=seq(-maxval,maxval,length=100),col=oceColorsFreesurface)

## training data
for(j in 1:4){
  png(paste0(file='plots/precip_global_',j,'.png'),width=1200,height=600)
  par(mgp = c(0,0,0), mar=c(0,0,0,0)) # bltr
  mapPlot(coastlineWorld, col="white", projection="+proj=moll",drawBox=FALSE)
  mapImage(unique(long),unique(lat),matrix(x[,j],nrow=288),colormap=cm)
  mapPolygon(coastlineWorld[['longitude']], coastlineWorld[['latitude']],
             border='grey5')
  dev.off()
}


## plot the legend
png(paste0(file='plots/precip_global_legend.png'),width=1200,height=600)
par(mgp = c(0,0,0), mar=c(0,0,0,0)) # bltr
drawPalette(colormap=cm)
dev.off()


## plot conditional samples
condsamp=samp.fixed; condsamp[ord,]=condsamp
for(j in 1:length(fi)) {
  png(paste0(file='plots/precip_global_cond_',fi[j],'.png'),width=1200,height=600)
  par(mgp = c(0,0,0), mar=c(0,0,0,0)) # bltr
  mapPlot(coastlineWorld, col="white", projection="+proj=moll",drawBox=FALSE)
  mapImage(unique(long),unique(lat),matrix(condsamp[,j],nrow=288),colormap=cm)
  mapPolygon(coastlineWorld[['longitude']], coastlineWorld[['latitude']],
             border='grey5')
  dev.off()
}

