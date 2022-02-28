# Scalable Bayesian transport maps for non-Gaussian spatial fields

To reproduce the figures and results in Katzfuss & Schäfer (2021), please follow the file `main.R`. Code for Figure 3 (written in Julia) can be found in the folder `figure3_julia`.

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

## fit posterior transport map
fit=optimFitMapAuto(Y,locs)

## draw new sample from posterior predictive distr
newsamp=condSampAuto(fit)

## plot new sample produced via transport map
quilt.plot(locs[,1],locs[,2],newsamp,nx=sqrt(N),ny=sqrt(N))
```

## Reference
Katzfuss, M. & Schäfer, F. (2021). Scalable Bayesian transport maps for high-dimensional non-Gaussian spatial fields. [*arXiv:2108.04211*](https://arxiv.org/abs/2108.04211).
