/*******************************************************************************
file pedometrics/src/MSSDCpp.cpp

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or 3 of the License
(at your option).

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.r-project.org/Licenses/

Purpose        : calculate the mean squared shortest distance used as criterion
                 to optimize spatial samples (spatial coverage sampling)
Author         : A. Samuel-Rosa <alessandrosamuelrosa at gmail.com>
Contributions  : 

Arguments:
x: matrix of distances between candidate locations and sample points
*******************************************************************************/
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".calcMSSDCpp")]]

double calcMSSDCpp(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow), xi(ncol);
  
  /* begin the main loop over the rows of the matrix of distances */
  for (int i = 0; i < nrow; i++) {
    
    /* begin secondary loop over the columns of the matrix of distances */
    /* every j-th element of the i-th row is copied to 'xi' */
    for (int j = 0; j < ncol; j++) {
      xi[j] = x(i, j);
    }

    /* get the minimum value in 'xi' and copy it to the i-th element of 'out'*/
    out[i] = min(xi);
  }
  
  /* square the minimim distances and return their mean */
  return mean(pow(out, 2));
}
/* End! */
