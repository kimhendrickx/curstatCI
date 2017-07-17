
<!-- README.md is generated from README.Rmd. Please edit that file -->
curstatCI
=========

The goal of curstatCI is to obtain confidence intervals for the distribution function of a random variable based on current status data. In the current status model, the variable of interest *X* with distribution function *F*<sub>0</sub> is not observed directly. A censoring variable *T* is observed instead together with the indicator *Δ* = (*X* ≤ *T*). *curstatCI* provides functions to estimate the distribution function *F*<sub>0</sub> and to construct pointswise confidence intervals around *F*<sub>0</sub>(*t*) based on an observed sample (*T*<sub>1</sub>, *Δ*<sub>1</sub>),…,(*T*<sub>*n*</sub>, *Δ*<sub>*n*</sub>) of size *n* from the observable random vector (*T*, *Δ*).

Installation
------------

You can install curstatCI from CRAN with:

``` r
# install.packages("curstatCI")
```

The package *curstatCI* requires the library *Rcpp*. To use the functions available in *curstatCI* load:

``` r
load(Rcpp)
load(curstatCI)
```

You can install curstatCI from github with:

``` r
# install.packages("devtools")
devtools::install_github("pietg/curstatCI")
```

Example
-------

This is a basic example which shows you how to obtain the confidence intervals for the distribution function of the time to infection for the Rubella data set. More information on the data and usage of the package can be found in the vignette "curstatCI"

``` r
library(Rcpp)
library(curstatCI)
set.seed(1)
data(rubella)
grid <-1:80
out<-ComputeConfIntervals(data=rubella,x=grid,alpha=0.05)
#> 
#> Piet Groeneboom and Kim Hendrickx 2017
#> For further information see:
#> The nonparametric bootstrap for the current status model,
#> Piet Groeneboom & Kim Hendricx, Electronical Journal of Statistics (2017).
#> 
#> The program produces the Studentized nonparametric bootstrap confidence intervals for the cdf, using the SMLE
#> 
#> 
#> Number of unique observations:    225
#> 
#> Sample size N =  230
#> 
#> 
#> The computations took    1.565   seconds

out$MLE
#>          [,1]      [,2]
#>  [1,]  0.0000 0.0000000
#>  [2,]  0.9452 0.2000000
#>  [3,]  1.2301 0.4857143
#>  [4,]  5.6411 0.5000000
#>  [5,]  9.4603 0.5714286
#>  [6,] 12.4548 0.8571429
#>  [7,] 15.4466 0.8666667
#>  [8,] 20.8000 0.8750000
#>  [9,] 23.2219 0.9285714
#> [10,] 25.2932 0.9401709
#> [11,] 77.8027 1.0000000

smle <-  out$SMLE
left<-out$CI[,1]
right<-out$CI[,2]

plot(grid,smle,type ='l', ylim=c(0,1), main= "",ylab="",xlab="",las=1)
lines(grid, left, col = 4)
lines(grid, right, col = 4)
segments(grid,left, grid, right)
```

<img src="README-unnamed-chunk-2-1.png" style="display: block; margin: auto;" />
