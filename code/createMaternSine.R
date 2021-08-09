### code that creates and samples from sine extension of matern


# ## provide the following settings:
# n=30^2
# Ns=c(10,20,50,100)
# N.test=50
# nonlin=2 # degree of nonlinearity
# reg=TRUE # regular grid or randomly sampled locations
# sepa=0 # separation for bimodal residual errors

m.max=30

## locations and ordering
grid.oneside=seq(0,1,length.out=sqrt(n))
if(reg) locs=as.matrix(expand.grid(grid.oneside,grid.oneside)) else
  locs=cbind(runif(n),runif(n))
ord=GPvecchia::order_maxmin_exact(locs)
locs.ord=locs[ord,]

## compute U and D from Matern
revMat=function(mat) mat[nrow(mat):1,ncol(mat):1,drop=FALSE]
range=.3; smooth=.5
cov.true=fields::Matern(rdist(locs.ord),range=range,nu=smooth)
U.unnormed=revMat(t(chol(solve(revMat(cov.true)))))
cond.sds=1/diag(U.unnormed)
U=t( t(U.unnormed)*cond.sds )

## coefficients for nearest neighbors
NNarray.max=GpGp::find_ordered_nn(locs.ord,m.max)[,-1]
weights=matrix(nrow=n,ncol=m.max)
for(i in 2:n) weights[i,]=-U[NNarray.max[i,],i]

## simulate from nonlinearized version
N.max=max(Ns)
N.total=N.max+N.test
x.total=array(dim=c(n,N.total))
mult=nonlin
f=function(vals,w){
  lin.pred=as.numeric(na.omit(vals*w))
  if(length(lin.pred)>2) temp=1:2 else if(length(lin.pred)>2) temp=c(1,1) else 
    temp=c()
  sum(lin.pred) + mult*sin(sum(lin.pred[temp])*4)
}
for(j in 1:N.total){
  x=numeric(length=n)
  for(i in 1:n){
    fx = f(x[NNarray.max[i,]],weights[i,])
    eps=rnorm(1,sample(sepa*c(-1,1),1)*cond.sds[i],cond.sds[i])
    x[i]=fx+eps
  }
  x.total[,j]=x
}
data.all=x.total[,1:N.max]
data.test=x.total[,(1+N.max):N.total]

score.exact=function(obs){
  ls=numeric(length=n)
  for(i in 1:n){
    fx = f(obs[NNarray.max[i,]],weights[i,])
    ls[i]=log.sum(log(.5)+dnorm(obs[i],fx+sepa*c(-1,1)*cond.sds[i],
                                cond.sds[i],log=TRUE))
  }
  sum(ls)
}

