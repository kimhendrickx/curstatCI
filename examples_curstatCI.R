################################################################
remove(list =ls())
library(Rcpp)
library(curstatCI)

# sample size
n <- 1000

# Exponential Data without ties (mean 1/lambda)
set.seed(100)
lambda = 0.75;
y <- rexp(n,lambda); summary(y)
t <- runif(n, 0,5)
delta <- as.numeric(y<=t)

# Normal Data N(20,5^2)
#set.seed(8455)
#y <- pmax(0, 20 + rnorm(n,0,5))
#t <- runif(n,0,40)
#delta <- as.numeric(y<=t)


## input data matrix
A<-cbind(t[order(t)], delta[order(t)], rep(1,n))


## input grid
maxt <- max(t)
maxgrid <-A[max(which(A[, 2]==0)),1]
mingrid <-A[min(which(A[, 2]==1)),1]
maxt;mingrid;maxgrid;



grid <-seq(mingrid,maxgrid ,length =200)
grid <-seq(A[2,1] ,A[n-1,1], length = 200)
grid <-seq(0,5, by = 0.05)

h =ComputeBW(A,grid)
plot(grid,h)

#smle <- ComputeSMLE(A, grid, h)
#smle <- ComputeSMLE(A, grid, rep(0.002, 200))
#plot(grid, smle, type = 'l')

out=ComputeConfIntervals(A,grid,0.05, h)
left = out$CI[,1]
right = out$CI[,2]


plot(grid, out$SMLE,type ='l')
lines(grid, left, col = 4)
lines(grid, right, col = 4)
segments(grid,left,grid,right)

#lines(grid,pnorm(grid,0.35,0.1), col =5)
#lines(grid,pnorm(grid,20,5), col =5)

lines(grid, 1-exp(-lambda*grid), col = 2)
abline(v = c(mingrid, maxgrid), col = 2 )



################################################################
## Rubella/HepatitisA data
remove(list =ls())
library(Rcpp)
library(curstatCI)

data(rubella)


grid <-seq(0.5, 79.5, by = 0.5)
#grid <- rubella$t
grid <- 1:85

set.seed(100)

#rubella
h=ComputeBW(rubella,grid)
plot(grid,h)

smle <- ComputeSMLE(rubella, grid, h)
plot(grid, smle, type = 'l')


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
grid <- 1:85

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

mingrid<-A[min(which(A[,2]>0)),1]
maxgrid<- A[max(which(A[,2]<A[,3])),1]
mingrid;maxgrid;

grid <-seq(0.02,1.98 ,by = 0.02)
ngrid <-length(grid)


h =ComputeBW(A,grid)
plot(grid,h)

out=ComputeConfIntervals(A,grid,0.05, h)
left = out$CI[,1]
right = out$CI[,2]


plot(grid, out$SMLE,type ='l', xlim = c(0, maxgrid))
lines(grid, left, col = 4)
lines(grid, right, col = 4)
lines(grid, (1-exp(-grid))/(1-exp(-2.0)), col = 2)


#####################

library(Rcpp)
library(curstatCI)

# sample size
n <- 1000

# Uniform data  U(0,2)
set.seed(2)
y <- runif(n,0,2)
t <- runif(n,0,2)
delta <- as.numeric(y <= t)

A<-cbind(t[order(t)], delta[order(t)], rep(1,n))
grid <-seq(0,2 ,by = 0.01)

#Bandwidth vector
h<-rep(2*n^-0.2,length(grid))


mle <- ComputeMLE(A)
smle <-ComputeSMLE(A,grid,h)
plot(grid, smle,type ='l', ylim=c(0,1), main= "",ylab="",xlab="",las=1)

