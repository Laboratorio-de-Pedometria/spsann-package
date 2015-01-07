/*******************************************************************************
file pedometrics/src/updatePPLCpp.cpp

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

Purpose        : Update the distance matrix used in the optimization of 
                 spatial samples using the number of points per lag distance 
                 class as criterion
Author         : A. Samuel-Rosa <alessandrosamuelrosa at gmail.com>
Contributions  : Joshua French <joshua.french at ucdenver.edu>

Arguments:
x: coordinates of all points
dm: matrix of distances
idx: row and column of the distance matrix to be updated
*******************************************************************************/
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".updatePPLCpp")]]

NumericMatrix updatePPLCpp(NumericMatrix x, NumericMatrix dm, int idx) {
  int ncol = x.ncol(), nrow = x.nrow(), i, j;
  NumericVector d(nrow, 0.0000);
  NumericMatrix x2(1, ncol);
  
  /* matrix indexes beging at 0 in C++, while in R it starts at 1 */
  idx -= 1;
  for (i = 0; i < ncol; i++) {
    x2(0, i) = x(idx, i);
  }
  
  /* begin the main loop over the rows of the matrix of coordinates */
  for (i = 0; i < nrow; i++) {
    
    /* begin the secondary loop over the columns of the matrix of coordinates */
    for (j = 0; j < ncol; j++) {
      d[i] += (x[nrow*j + i] - x2[j])*(x[nrow*j + i] - x2[j]);
    }
    
    /* take the squared root */
    d[i] = pow(d[i], 0.5);
  }
  
  /* replace the values in the distance matrix */
  for (i = 0; i < nrow; i++) {
    dm(i, idx) = d[i];
    dm(idx, i) = d[i];
  }
  return (dm);
}
/* End! */
