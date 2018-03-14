# Releasing package curstatCI
12-10-2017

## Test environments 

* Mac environment  (macOS Sierra version 10.12.5)
    * platform: x86_64-apple-darwin15.6.0 (64-bit) 
    * R version 3.4.1  (2017-06-30) 


* Windows environment  (Microsoft Windows 10 Enterprise)
    * platform: x86_64-w64-mingw32/x64 (64-bit) 
    * R version 3.4.1 Patched (2017-07-01 r72876)


* Linux environment (Ubuntu 16.04 LTS)
    * platform: x86_64-pc-linux-gnu (64-bit) 
    * R version 3.4.2 (2017-09-28) 

    
* win-builder
    * platform: x86_64-w64-mingw32 (64-bit)
    * R Under development (unstable) (2017-09-12 r73242)

## R CMD check results    

There are no errors, no warnings and no notes.


## Reverse Dependencies
There are no problems with reverse dependencies.


## Resubmission
This is a resubmission. In this version I have: 

* corrected an error in the source code of the function ComputeBW.

* changed the reference "Groeneboom and Hendrickx (2017) <arXiv:1701.07359>" 
    * to "Groeneboom and Hendrickx (2017) <doi:10.1214/17-EJS1345"" in the DESCRIPTION 
    * to "Groeneboom, P. and Hendrickx, K. (2017). The nonparametric bootstrap for the current status model. Electronic Journal of Statistics 11(2):3446-3848." in the references of the documentation of the functions ComputeMLE, ComputeSMLE, ComputeBW and ComputeConfIntervals and the references of the vignette curstatCI. 
