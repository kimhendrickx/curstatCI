#include "curstat.h"

//////////////////////////////////////////////////////////////////////////////////
//
// Confidence Intervals
//
// Export: Yes
// Date Created: 11.07.2017
// Created by: Piet Groeneboom and Kim Hendrickx
//////////////////////////////////////////////////////////////////////////////////
//'@title Pointwise Confidence Intervals under Current Status data
//'
//'@description The function ComputeConfIntervals computes pointwise confidence intervals for the distribution function under current status data.
//'The confidence intervals are based on the Smoothed Maximum likelihood Estimator and constructed using the nonparametric bootstrap.
//'
//'@param data Dataframe with three variables:
//'\describe{
//'     \item{t}{Observation points t sorted in ascending order. All observations need to be positive. The total number of unique observation points equals \code{length(t)}.}
//'     \item{freq1}{Frequency of observation t satisfying \eqn{x \le t}.
//'                   The total number of observations with censoring indicator \eqn{\delta =1} equals \code{sum(freq1)}. }
//'     \item{freq2}{Frequency of observation t. The sample size equals \code{sum(freq2)}. If no tied observations are present in the data \code{length(t)} equals \code{sum(freq2)}. }
//'}
//'
//'@param x numeric vector
//'
//'@param alpha confidence level of pointwise confidence intervals
//'
//'@return List with 4 variables:
//'
//'\enumerate{
//'     \item{MLE }{ Maximum Likelihood Estimator}
//'     \item{SMLE }{ Smoothed Maximum Likelihood Estimator in x }
//'     \item{h }{ data-driven bandwidth choice for each point in x}
//'     \item{CI }{ pointwise confidence interval for each points in x}
//' }
//'
//' @details In the current status model, the variable of interest \eqn{X} with distribution function \eqn{F} is not observed directly.
//'A censoring variable \eqn{T} is observed instead together with the indicator \eqn{\Delta = (X \le T)}.
//' ComputeConfIntervals computes the pointwise \code{1-alpha} bootstrap confidence intervals around the SMLE of \eqn{F} based on a sample of size \code{n <- sum(data$freq2)} from the observable random  vector \eqn{(T, \Delta)}.
//' The bandwidth used for estimating the SMLE at a point in the vector x is obtained by minimizing the pointwise Mean Squared Error using the subsampling pricinciple in combination with undersmoothing.
//'
//' @seealso \code{vignette("curstatCI")}
//'
//'@examples
//'library(Rcpp)
//'library(curstatCI)
//'
//'# sample size
//'n <- 1000
//'
//'# Uniform data  U(0,2)
//'set.seed(2)
//'y <- runif(n,0,2)
//'t <- runif(n,0,2)
//'delta <- as.numeric(y <= t)
//'
//'A<-cbind(t[order(t)], delta[order(t)], rep(1,n))
//'grid<-seq(0,2 ,by = 0.01)
//'
//'# pointwise confidence intervals at grid points:
//'out<-ComputeConfIntervals(data = A,x =grid,alpha = 0.05)
//'
//'left <- out$CI[,1]
//'right <- out$CI[,2]
//'
//'plot(grid, out$SMLE,type ='l', ylim=c(0,1), main= "",ylab="",xlab="",las=1)
//'lines(grid, left, col = 4)
//'lines(grid, right, col = 4)
//'segments(grid,left, grid, right)
//'
//'@references The nonparametric bootstrap for the current status model, Groeneboom, P. and Hendrickx, K. Electronical Journal of Statistics (2017)
//'@export
// [[Rcpp::export]]
List ComputeConfIntervals(DataFrame data, NumericVector x, double alpha)
{
  double          A,B,*data0,*Data,*data1,*data2,*grid,*p,*p2,h;
  double          *tt,**f3,*lowbound,*upbound,*f4, **MSE, *hmin1;
  double          *cumw,*cs,*F,*F2,*jumploc,*y,*y2,*SMLEa ,*SMLE,*SMLE2a, *SMLE2,*cc;
  int             i,j,k,m,m2,N,n,*delta,**freq,njumps,*delta2,*freq1,*freq2,**frequence1;
  int             percentile1,percentile2,iter,iter2,NumIt=1000,nB, npoints, nIterH;
  clock_t         StartTime, StopTime;
  double          Time_bootstrap;


  DataFrame DF = Rcpp::DataFrame(data);
  NumericVector xx = DF[0];
  IntegerVector freq01 = DF[1];
  IntegerVector freq02 = DF[2];

  Rcout << std::endl;
  //Rcout << "Piet Groeneboom and Kim Hendrickx 2017" << std::endl << "For further information see:" << std::endl;
  Rcout << "The program produces the Studentized nonparametric bootstrap confidence intervals for the cdf, using the SMLE" << std::endl<< std::endl;
  Rcout << "For further information see:" << std::endl;
  Rcout << "The nonparametric bootstrap for the current status model," << std::endl;
  Rcout << "Piet Groeneboom & Kim Hendrickx, Electronical Journal of Statistics (2017)." << std::endl << std::endl;


  // determine the number of rows of the data frame

  percentile1=(int)((alpha/2)*NumIt);
  percentile2=(int)((1-alpha/2)*NumIt);

  n = (int)xx.size();
  npoints = (int)x.size();


  data0= new double[n+1];
  freq1= new int[n+1];
  freq2= new int[n+1];

  data0[0]=0;

  for (i=1;i<=n;i++)
  {
    data0[i]=(double)xx[i-1];
    freq1[i]=(int)freq01[i-1];
    freq2[i]=(int)freq02[i-1];
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  A = 0.0;
  B = data0[n];



  grid = new double[npoints+1];
  grid[0] = A;

  for (i=1;i<=npoints;i++)
    grid[i]= x[i-1];
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  N=0;
  for (i=1;i<=n;i++)
    N += freq2[i];


  nB = (int) N/2;
  if(N>= 100 && N < 1000)
    nB = 50;
  if(N>= 1000 && N < 10000)
    nB = 200;
  if(N>= 10000)
    nB = 500;

  Rcout << "Number of unique observations:" << std::setw(7) << n << std::endl;
  Rcout << "Sample size n =" << std::setw(5) << N << std::endl;
  //Rcout << "Left point of observation interval A = " <<  std::setw(5) << A << std::endl;
  //Rcout << "Right point of observation interval B = " <<  std::setw(5) << B << std::endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // form the data vector for the bootstrap

  Data = new double[N+1];
  delta = new int[N+1];
  tt=new double[N+1];
  delta2 = new int[N+1];

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // grid search for optimal bandwidth

  nIterH = 100;
  cc = new double[nIterH+1];
  for (i=0;i<=nIterH;i++)
    cc[i]= 2.5*(B-A)*i/nIterH;

  MSE = new double*[npoints+1];
  for (iter2=0;iter2<=npoints;iter2++)
    MSE[iter2] = new double[nIterH+1];

  hmin1 = new double[npoints+1];

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  Data[0]=0;

  j=0;

  for (i=1;i<=n;i++)
  {
    for (k=1;k<=freq1[i];k++)
    {
      j++;
      Data[j]=data0[i];
      delta[j]=1;

    }
    for (k=1;k<=freq2[i]-freq1[i];k++)
    {
      j++;
      Data[j]=data0[i];
      delta[j]=0;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  F= new double[n+1];
  F2= new double[n+1];
  cumw= new double[n+1];
  cs= new double[n+1];
  y= new double[n+1];
  y2= new double[n+1];
  jumploc= new double[n+1];
  data1= new double[n+1];
  data2= new double[n+1];

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  SMLE= new double[npoints+1];
  SMLE2= new double[npoints+1];
  SMLEa= new double[npoints+1];
  SMLE2a= new double[npoints+1];
  p= new double[n+1];
  p2= new double[n+1];

  frequence1 = new int*[2];
  for (i=0;i<2;i++)
    frequence1[i] = new int[n+1];

  freq = new int*[2];
  for (i=0;i<2;i++)
    freq[i] = new int[n+1];

  f3  = new double*[NumIt+1];

  for (iter=0;iter<NumIt+1;iter++)
    f3[iter] = new double[npoints];


  F[0]=F2[0]=0;
  cumw[0]=cs[0]=0;

  y[0]=y2[0]=0;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  StartTime = clock();

  for (i=1;i<=n;i++)
  {
    cs[i]=cs[i-1]+(double)freq1[i];
    cumw[i]=cumw[i-1]+(double)freq2[i];
    frequence1[1][i]=freq1[i];
    frequence1[0][i]=freq2[i]-freq1[i];
  }


  convexmin(n,cumw,cs,y);

  j=0;

  jumploc[0]=0;

  for (i=1;i<=n;i++)
  {
    if (y[i]>y[i-1])
    {
      j++;
      p[j]=y[i]-y[i-1];
      F[j]=y[i];
      jumploc[j]=data0[i];
    }
  }

  njumps=j;

  NumericMatrix out1 = NumericMatrix(njumps+1,2);

  for (i=0;i<=njumps;i++)
  {
    out1(i,0)=jumploc[i];
    out1(i,1) = F[i];
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////



  // bandwidth for initial SMLE
  h = 2*(B-A)*pow(N,-1.0/5);

  SMLEa[0]=0;

  for (i=1;i<=npoints;i++)
    SMLEa[i]=bdf(A,B,njumps,jumploc,p,grid[i], h);


  for (i=1;i<=npoints;i++)
  {
    for (j=1;j<=nIterH;j++)
    {
      MSE[i][j]=0;
    }
  }




  ////////////////////////////////////////////////////////////////////////////////////////////////////

  for (iter=0;iter<NumIt;iter++)
  {
    data_bootstrap2(N,nB,n,&m,Data,tt,data1,freq,delta,delta2);

    for (i=1;i<=m;i++)
    {
      cs[i]=cs[i-1]+(double)freq[1][i];
      cumw[i]=cumw[i-1]+(double)(freq[0][i]+freq[1][i]);
    }

    convexmin(m,cumw,cs,y2);


    j=0;

    for (i=1;i<=m;i++)
    {
      if (y2[i]>y2[i-1])
      {
        j++;
        data2[j]=data1[i];
        p2[j]=y2[i]-y2[i-1];
        F2[j]=y2[i];
      }
    }

    m2=j;


    SMLE2a[0]=0;


    for (i=1;i<=npoints;i++)
    {
      for (j=1;j<=nIterH;j++)
      {
        SMLE2a[i]=bdf(A,B,m2,data2,p2,grid[i],cc[j]*pow(nB,-1.0/5));

        MSE[i][j] += SQR(SMLE2a[i]-SMLEa[i])/NumIt;
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Select the optimal bandwidth and bias estimate

  for (i=1;i<=npoints;i++)
  {
    hmin1[i]=cc[1]*pow(N,-1.0/4);

    for (j=2;j<=nIterH;j++)
    {
      if (MSE[i][j]<MSE[i][j-1])
        hmin1[i]=cc[j]*pow(N,-1.0/4);
    }
  }


  // ## Compute SMLE  with optimal bandwidth
  SMLE[0]=0;

  for (i=1;i<=npoints;i++)
  {
    SMLE[i]=bdf(A,B,njumps,jumploc,p,grid[i],hmin1[i]);
  }

  NumericVector out2(npoints);
  NumericVector out3(npoints);

  for (i=0;i<npoints;i++)
  {
    out2[i] = SMLE[i+1];
    out3[i] = hmin1[i+1];
  }


  for (iter=0;iter<NumIt;iter++)
  {
    data_bootstrap(N,n,&m,Data,tt,data1,freq,delta,delta2);

    for (i=1;i<=m;i++)
    {
      cs[i]=cs[i-1]+(double)freq[1][i];
      cumw[i]=cumw[i-1]+(double)(freq[0][i]+freq[1][i]);
    }

    convexmin(m,cumw,cs,y2);


    j=0;

    for (i=1;i<=m;i++)
    {
      if (y2[i]>y2[i-1])
      {
        j++;
        data2[j]=data1[i];
        p2[j]=y2[i]-y2[i-1];
        F2[j]=y2[i];
      }
    }

    m2=j;


    SMLE2[0]=0;


    for(i =1; i<=npoints; i++)
    {
      SMLE2[i]=bdf(A,B,m2,data2,p2,grid[i],hmin1[i]);

      if(sqrt(varF(N,n,frequence1,y,data0[1],data0[n],data0,hmin1[i],grid[i])) > 0)
      {
        f3[iter][i]=(SMLE2[i]-SMLE[i])/( (sqrt( varF(N,m,freq,y2,0.0,B,data1,hmin1[i],grid[i]))));
      }
      else{
        f3[iter][i]=(SMLE2[i]-SMLE[i])/( (sqrt( varF(N,m,freq,y2,0.0,B,data1,(B-A)/2,grid[i]))));
      }

    }

    //Rcout << "Bootstrap iteration " << setprecision(6) << iter +1 << " out of " << setprecision(6)  << NumIt << std::endl;

  }


  StopTime  = clock();
  Time_bootstrap   = (double)(StopTime - StartTime)/(double)CLOCKS_PER_SEC;

  f4= new double[NumIt+1];

  lowbound=new double[npoints];
  upbound=new double[npoints];

  for (i=1;i<=npoints;i++)
  {
    for (iter=0;iter<NumIt;iter++)
      f4[iter]=f3[iter][i];

    qsort(f4,NumIt,sizeof(double),compare);

    lowbound[i]= fmax(0,SMLE[i]-f4[percentile2-1]*sqrt(varF(N,n,frequence1,y,data0[1],data0[n],data0,hmin1[i],grid[i])));;
    upbound[i]= fmin(1,SMLE[i]-f4[percentile1-1]*sqrt(varF(N,n,frequence1,y,data0[1],data0[n],data0,hmin1[i],grid[i])));;
  }

  if(grid[1] == 0)
    upbound[1] = 0;

  Rcout << std::endl << std::endl;
  Rcout << "The computations took    " << setprecision(10) << Time_bootstrap << "   seconds"  << std::endl;

  NumericMatrix out4 = NumericMatrix(npoints,2);

  for (i=0;i<npoints;i++)
  {
    out4(i,0)=lowbound[i+1];
    out4(i,1)=upbound[i+1];
  }

  Rcout << std::endl;

  /*
  ofstream file_("CI_SMLE.txt");

  if (file_.is_open())
  {
  for (i=1;i<npoints;i++)
  {
  file_ << setprecision(10) << setw(20) << grid[i];
  file_ << setprecision(11) <<  setw(20) << lowbound[i]
        << setprecision(11) <<  setw(20) << upbound[i];
  file_ << "\n";
  }
  file_.close();
  }
  */



  // make the list for the output, containing the MLE, hazard, the bootstrap confidence intervals and -log likelihood

  List out = List::create(Rcpp::Named("MLE")=out1,Rcpp::Named("SMLE")=out2,Rcpp::Named("h")=out3,Rcpp::Named("CI")=out4 );

  // free memory

  delete[] Data, delete[] delta, delete[] tt, delete[] delta2, delete[] SMLE, delete[] SMLE2,
  delete[] F, delete[] F2, delete[] cumw, delete[] cs, delete[] y, delete[] y2,
  delete[] jumploc,  delete[] data0, delete[] data1, delete[] data2, delete[] p, delete[] p2,
  delete[] lowbound, delete[] upbound;

  for (i = 0;i<2;i++)
    delete[] freq[i];
  delete[] freq;

  for (iter = 0;iter < NumIt;iter++)
    delete[] f3[iter];
  delete[] f3;

  return out;
}
