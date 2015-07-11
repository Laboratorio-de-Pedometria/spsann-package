## Test environments
* local x86_64-pc-linux-gnu (ubuntu 14.), R 3.2.1
* ubuntu 12.04 (on travis-ci), R 3.2.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE in win-builder (devel).

   Undefined global functions or variables:
      axis cor hist legend lines mtext par plot points quantile runif
      setTxtProgressBar terms text txtProgressBar

This note is due to changes in R-devel utilities on 29 Jun 2015. According to 
R-devel/NEWS, `R CMD check --as-cran` now checks code usage (via `codetools`) 
with only the `base` package attached, so that functions from default packages 
other than `base` which are used in the package code but not imported are 
reported as undefined globals.
Source: http://developer.r-project.org/blosxom.cgi/R-devel/2015/06/29#n2015-06-29

## Downstream dependencies
There are no downstream dependencies of `spsann`
