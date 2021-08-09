#######   R code for Katzfuss & Schaefer  #########

# This file contains or sources the code to reproduce all plots and results.
# Set your working directory to the folder containing this file.
# Place all additional R files in a subfolder called "code/".
# Place prec_days.RData in a subfolder called "data/".
# Ensure that there are (empty) subfolders called "output/" and "plots/".
# Note: Some portions of the code can take a very long time to run.


### required functions
# setwd("G:/My Drive/projects/triangular_spatial")
source('code/nonlinearSpatial_functions.R')
source('code/logScore_comparison.R')



######  illustration of maximin ordering  ######

## compute ordering, NN, and length scales
n=60^2
Ns=100; N.test=1
nonlin=0 # degree of nonlinearity
reg=TRUE # regular grid or randomly sampled locations
sepa=0 # separation for bimodal residual errors
source('code/createMaternSine.R')
NNarray=NNarray.max[,1:4]
scales=computeScales(locs.ord,NNarray)

## plot orderings
is=c(13,51,290)
for(i in is){
  pdf(file=paste0('plots/maxmin',i,'.pdf'),width=4.0,height=4.0)
  par(mgp = c(1.6,.5,0), mar=c(.1,.1,.1,.1)) # bltr
  plot(locs.ord[,1],locs.ord[,2],pch='.',xaxt='n',yaxt='n',
       xlab='',ylab='',col='grey',cex=1.5)
  points(locs.ord[1:(i-1),1],locs.ord[1:(i-1),2],col=1,cex=1.5)
  points(locs.ord[NNarray[i,],1],locs.ord[NNarray[i,],2],
         col=3,pch='x',cex=2)
  points(locs.ord[i,1],locs.ord[i,2]-1/80,col=4,pch='+',cex=3)
  segments(locs.ord[i,1],locs.ord[i,2],locs.ord[NNarray[i,1],1],
           locs.ord[NNarray[i,1],2],col=2,lwd=2)
  if(i==13) text(mean(locs.ord[c(i,NNarray[i,1]),1])+.03,
                 mean(locs.ord[c(i,NNarray[i,1]),2])-.03,
                 expression(italic(l[i])),cex=2,col=2)
  dev.off()
}

## plot length scales (min distances)
pdf(file='plots/scales.pdf',width=3.0,height=3.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
plot(1:n,scales,xlab='i',ylab='min dist',cex=1,col=2)
lines((1:n)^(-.5),col=4)
legend('topright',c(expression(min~dist~~italic(l[i])),expression(1/sqrt(i))),
       pch=c(1,NA),lty=c(NA,1),col=c(2,4),pt.cex=c(1,NA))
dev.off()



### decay of Matern coefficients

# conditional sd
theta=coef(lm(log(cond.sds) ~ log(scales)))

pdf(file='plots/matern_nugget.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
plot(cond.sds,xlab='i',ylab='residual SD')
lines(exp(theta[1]+theta[2]*log(scales)),col=2)
dev.off()


# plot of kriging weights
pdf(file='plots/matern_sq_weights.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
squared.dev=apply(weights^2,2,mean,na.rm=TRUE)
plot(squared.dev[1:15],xlab='k',ylab='mean squared weights',col=1,cex=1.5)
abline(h=0)
lines(squared.dev[1:15],col=1,lwd=2)
dev.off()



### scatterplot vs NN
i=290
k=5
NN=NNarray.max[i,]
ra=range(data.all[NN[c(1,k)],])

pdf(file='plots/matern3d.pdf',width=4.0,height=4.0)
par(mgp = c(.1,.1,0), mar=c(.1,.1,.1,.1)) # bltr
s3d=scatterplot3d::scatterplot3d(t(data.all[c(NN[1],NN[k],i),]),angle=30,
                                 color=rainbow(100)[cut(c(data.all[i,]),100)],
                                 xlim=ra,ylim=ra,pch=20,xlab='1st NN',
                                 zlab=expression(y[i]),ylab='', #paste0(k,'th NN'),
                                 cex.symbols = 1.2)
                             # type='h',box=FALSE,highlight.3d=TRUE)
s3d$box3d(col='grey')
s3d$plane3d(c(0, #sum(weights[i,-c(1,k)]*rowMeans(data.all[NN[-c(1,k)],])),
              weights[i,c(1,k)]),col=1)
dims <- par("usr")
x <- dims[1]+ 0.85*diff(dims[1:2])
y <- dims[3]+ 0.08*diff(dims[3:4])
text(x,y,paste0(k,'th NN'),srt=30)
dev.off()








######  matern+sine illustrations  ######


source('code/maternSine_illustration.R')







######  matern+sine KL comparisons  ######

name=c('linear','S-linear','nonlin','S-nonlin','DPM')

### compute KL values for different settings, save in output folder
source('code/maternSine_comparison.R')

### plot comparisons results for increasing N
files=c('sin_lin','sin','sin_random','sin_DPM')
for(fi in 1:length(files)){
  load(file=paste0('output/compRes_',files[fi],'.RData'))
  pdf(file=paste0('plots/klcomp_',files[fi],'.pdf'),width=3.5,height=3.5)
  par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
  matplot(log10(Ns),kls,type='l',lwd=3.5,xlab='n',ylab='KL',xaxt='n',
          ylim=range(kls[,-4],na.rm=TRUE)) # exclude nonlin for plot range
  axis(1,at=log10(Ns),labels=Ns,lwd=0,lwd.ticks=1)
  if(fi %in% c(1,4)) lloc='topright' else lloc='bottomleft' 
  legend(lloc,name[methods],lty=methods,col=methods,lwd=2.5,bg='white')
  dev.off()
}




######  prepare precip anomalies for each day in june  ######

load(file='data/prec_days.RData') # precs=[locs,98yrs,30d]
locs=cbind(lon,lat)

## log transform
logprecs=log(precs+1e-10)

## anomalies: standardize according to june 1
mus=apply(logprecs[,,1],1,mean)
sds=apply(logprecs[,,1],1,sd)
xs=array(dim=dim(logprecs))
for(d in 1:30) xs[,,d]=t(scale(t(logprecs[,,d]),mus,sds))

save(xs,locs,file='output/prec_anomalies.RData')




######  climate data illustrations  ######

source('code/climate_illustrations.R')




######  climate data: log-score comparison  ######

source('code/climate_comparison.R')

load(file='output/compRes_climate_scales.RData')
als=-apply(ls,2:3,mean,na.rm=TRUE)

pdf(file='plots/ls_precip.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
matplot(Ns,als,type='l',lwd=2,xlab='n',ylab='LS',
        ylim=range(als[,c(1,3,5)]),col=c(1:6,8))
legend('topright',c(name,'tapSamp','expCov'),lty=1:7,
       col=c(1:6,8),bg='white',lwd=2)
dev.off()
round(apply(als,2,min))



######  comparison for global climate data  ######

source('code/climate_comparison_global.R')

load(file='output/compRes_climate_global.RData')
als=-apply(ls,2:3,mean,na.rm=TRUE)

pdf(file='plots/ls_precip_global.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
matplot(Ns,als,type='l',lwd=2,xlab='n',ylab='LS',
        lty=methods,col=methods)
legend('topright',name[methods],lty=methods,col=methods,lwd=2.5)
dev.off()



######  plot global data and conditional samples  ######

source('code/draw_global_maps.R')



######  examine fit for largest setting  ######

source('code/climate_saveTheta.R')

load(file='output/fit_global.RData')

n
m.threshold(fg$theta,m.max=30)

## posterior median of d_i
d.med=sqrt(qinvgamma(.5,fg$alpha.posts,fg$beta.posts))
sum(d.med>.05*max(d.med))
# plot(d.med,xlab='maximin index',
#      ylab=expression(posterior~median~of~d[i]))
# plot(sort(d.med,TRUE))
