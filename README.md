# Scalable Bayesian transport maps for non-Gaussian spatial fields

To reproduce the figures and results in Katzfuss & Schäfer (2021), please follow the file `main.R`.

To get started using the method on your own data, the following toy example may be helpful.
```{r}
## load the "package"
source('https://raw.githubusercontent.com/katzfuss-group/BaTraMaSpa/main/code/nonlinearSpatial_functions.R')

## simulate some toy data
N=20^2; n=50
grid=seq(0,1,length.out=sqrt(N))
locs=as.matrix(expand.grid(grid,grid))
cov.chol=t(chol(exp(-rdist(locs)/.3)))
Y=cov.chol%*%matrix(rnorm(N*n),nrow=N)

## ordering and nearest neighbors
ord=GPvecchia::order_maxmin_exact(locs)
NNarray.max=GpGp::find_ordered_nn(locs[ord,],30)[,-1]
scales=computeScales(locs[ord,],NNarray.max)
Y.ord=Y[ord,]

## fit posterior transport map
fit=optimFitMap(Y.ord,NNarray.max,scales=scales)

## draw new sample from posterior predictive distr
newsamp=condSamp(fit)
newsamp[ord]=newsamp  # original ordering

## plot the new sample
quilt.plot(locs[,1],locs[,2],newsamp,nx=sqrt(N),ny=sqrt(N))
```

## Reference
Katzfuss, M. & Schäfer, F. (2021). Scalable Bayesian transport maps for high-dimensional non-Gaussian spatial fields. [*arXiv:2108.04211*](https://arxiv.org/abs/2108.04211).
