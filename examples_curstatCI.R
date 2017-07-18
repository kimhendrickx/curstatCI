################################################################
remove(list =ls())
library(Rcpp)
library(curstatCI)


n <- 1000



# Exponential Data without ties (mean 1/lambda)
set.seed(100)
lambda = 2;
y <- rexp(n,lambda)
t <- rexp(n,lambda)

# Normal Data N(0.35,0.01)
set.seed(100)
y <- pmax(0, 0.35+ rnorm(n,0,0.1))
t <- runif(n,0,2)


delta <- as.numeric(y<=t)

A<-cbind(t[order(t)], delta[order(t)], rep(1,n))
grid <-seq(0.02,1.98 ,by = 0.02)
ngrid <-length(grid)


h =ComputeBW(A,grid)
plot(grid,h)

out=ComputeConfIntervals(A,grid,0.05, h)
left = out$CI[,1]
right = out$CI[,2]

plot(grid, out$SMLE,type ='l')
lines(grid, left, col = 4)
lines(grid, right, col = 4)
segments(grid,left,grid,right)




################################################################
## Rubella/HepatitisA data
remove(list =ls())
library(Rcpp)
library(curstatCI)

data(rubella)


grid <-1:75
#grid <- rubella$t

set.seed(100)

#rubella
h=ComputeBW(rubella,grid)
plot(grid,h)

#out=ComputeConfIntervals(rubella,grid,0.05, rep(20,85))
out=ComputeConfIntervals(rubella,grid,0.05, h)

left = out$CI[,1]
right = out$CI[,2]

plot(grid, out$SMLE,type ='l')
lines(grid, left, col = 4)
lines(grid, right, col = 4)
segments(grid,left,grid,right)

#hepatitisA
data(hepatitisA)


h=ComputeBW(hepatitisA,grid)
plot(grid,h)

out=ComputeConfIntervals(hepatitisA,grid,0.05, h)
left = out$CI[,1]
right = out$CI[,2]

plot(grid, out$SMLE,type ='l')
lines(grid, left, col = 4)
lines(grid, right, col = 4)
segments(grid,left,grid,right)

################################################################
#Truncated exponential
remove(list =ls())
library(Rcpp)
library(curstatCI)
set.seed(100)

n<-1000
t <- rep(NA, n)
delta <- rep(NA, n)

for(i in (1:n) ){
  x= runif(1);
  y=-log(1-(1-exp(-2))*x);
  t[i]=2*runif(1);
  if(y<=t[i]){
    delta[i]=1
  }
  else{
    delta[i]=0
  }
}

A=cbind(t[order(t)], delta[order(t)], rep(1,n))
delta <- as.numeric(y<=t)

grid <-seq(0.02,1.98 ,by = 0.02)
ngrid <-length(grid)


h =ComputeBW(A,grid)
plot(grid,h)

out=ComputeConfIntervals(A,grid,0.05, h)
left = out$CI[,1]
right = out$CI[,2]


plot(grid, out$SMLE,type ='l')
lines(grid, left, col = 4)
lines(grid, right, col = 4)
lines(grid, (1-exp(-grid))/(1-exp(-2.0)), col = 2)
