## ------------------------------------------------------------------------
set.seed(100)

n<-1000
t<-rep(NA, n)
delta<-rep(NA, n)

for(i in (1:n) ){
  x<-runif(1)
  y<--log(1-(1-exp(-2))*x)
  t[i]<-2*runif(1);
  if(y<=t[i]){ delta[i]<-1}
  else{delta[i]<-0}}

## ------------------------------------------------------------------------
A<-cbind(t[order(t)], delta[order(t)], rep(1,n))
head(A)

## ---- fig.width=4,fig.height=4, fig.align='center'-----------------------
library(Rcpp)
library(curstatCI)

mle<-ComputeMLE(data=A)
mle

## ---- fig.width=4,fig.height=4, fig.align='center'-----------------------
plot(mle$x, mle$mle,type='s', ylim=c(0,1),xlim=c(0,2), main="",ylab="",xlab="",las=1)

## ---- fig.width=4,fig.height=4, fig.align='center'-----------------------
grid<-seq(0.02,1.98, 0.02)
bw<-rep(2*n^-0.2,length(grid))
smle<-ComputeSMLE(data=A, x=grid, bw=bw)
plot(grid, smle,type='l', ylim=c(0,1), main="",ylab="",xlab="",las=1)
lines(grid, (1-exp(-grid))/(1-exp(-2.0)), col=2, lty=2)

## ------------------------------------------------------------------------
c(min(t),max(t))
grid<-seq(0.01,1.99 ,by=0.01)

## ---- fig.width=4,fig.height=4, fig.align='center'-----------------------
bw<-ComputeBW(data=A, x=grid)
plot(grid, bw, main="",ylim=c(0.5,0.7),ylab="",xlab="",las=1)

## ---- fig.width=4,fig.height=4, fig.align='center'-----------------------
out<-ComputeConfIntervals(data=A,x=grid,alpha=0.05, bw=bw )

## ------------------------------------------------------------------------
attributes(out)

out$NonStudentized

## ---- fig.width=4,fig.height=4, fig.align='center'-----------------------
left<-out$CI[,1]
right<-out$CI[,2] 

plot(grid, out$SMLE,type='l', ylim=c(0,1), main="",ylab="",xlab="",las=1)
lines(grid, left, col=4)
lines(grid, right, col=4)
segments(grid,left, grid, right)
lines(grid, (1-exp(-grid))/(1-exp(-2.0)), col=2)

## ------------------------------------------------------------------------
data(hepatitisA)
head(hepatitisA)

## ---- results='hide'-----------------------------------------------------
grid<-1:75

bw<-ComputeBW(data=hepatitisA,x=grid)
out<-ComputeConfIntervals(data=hepatitisA,x=grid,alpha=0.05, bw=bw)

## ---- fig.width=4,fig.height=4, fig.align='center'-----------------------
out$SMLE[18]
left<-out$CI[,1]
right<-out$CI[,2]

plot(grid, out$SMLE,type='l', ylim=c(0,1), main="",ylab="",xlab="",las=1)
lines(grid, left, col=4)
lines(grid, right, col=4)
segments(grid,left, grid, right)

## ------------------------------------------------------------------------
out$NonStudentized

## ------------------------------------------------------------------------
data(rubella)
head(rubella)
summary(rubella$t)

## ---- results="hide"-----------------------------------------------------
grid<-1:80
bw<-ComputeBW(data=rubella,x=grid)
out<-ComputeConfIntervals(data=rubella,x=grid,alpha=0.05, bw=bw)

## ---- fig.width=4,fig.height=4, fig.align='center'-----------------------
left<-out$CI[,1]
right<-out$CI[,2]

plot(grid, out$SMLE,type='l', ylim=c(0,1), main="",ylab="",xlab="",las=1)
lines(grid, left, col=4)
lines(grid, right, col=4)
segments(grid,left, grid, right)

## ---- fig.width=4,fig.height=4, fig.align='center'-----------------------
# sample size
n<-1000

# Exponential Data without ties (mean 1/lambda)
set.seed(100)
lambda<-0.75;
y<-rexp(n,lambda); summary(y)
t<-runif(n, 0,5)
delta<-as.numeric(y<=t)

# input data matrix
A<-cbind(t[order(t)], delta[order(t)], rep(1,n))


# x vector (recommended)
maxgrid<-A[max(which(A[, 2]==0)),1]
mingrid<-A[min(which(A[, 2]==1)),1]
c(mingrid,maxgrid)

# x vector 2
grid <-seq(0,5, by = 0.05)

bw<-ComputeBW(A,grid)
plot(grid,bw,ylim=c(1,2), main="",ylab="",xlab="",las=1)


out=ComputeConfIntervals(A,grid,0.05, bw)
left<-out$CI[,1]
right<-out$CI[,2]

## ---- fig.width=4,fig.height=4, fig.align='center'-----------------------
plot(grid, out$SMLE,type='l', ylim=c(0,1), main="",ylab="",xlab="",las=1)
lines(grid, left, col=4)
lines(grid, right, col=4)
segments(grid,left,grid,right)

# True distribution function
lines(grid, 1-exp(-lambda*grid), col=2)

# Recommended x limits
abline(v=c(mingrid, maxgrid), col=2 )

