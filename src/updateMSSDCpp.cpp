/*******************************************************************************
file pedometrics/src/updateMSSDCpp.cpp

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

Purpose        : Update the matrix of distances used in the optimization of 
                 spatial samples using the mean squared shortest distance as 
                 criterion
Author         : A. Samuel-Rosa <alessandrosamuelrosa at gmail.com>
Contributions  : Joshua French <joshua.french at ucdenver.edu>

Arguments:
x1: coordinates of the prediction points
x2: coordinates of the new point
dm: matrix of distances
idx: column of the matrix of distances to be updated
*******************************************************************************/
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".updateMSSDCpp")]]

NumericMatrix updateMSSDCpp(NumericMatrix x1, NumericMatrix x2, 
                            NumericMatrix dm, int idx) {
  int ncol = x1.ncol(), nrow = x1.nrow(), i, j;
  
  /* Fill the vector with zeros and four decimal places to set the format of */
  /* the output */
  NumericVector d(nrow, 0.0000);
  
  /* begin the main loop over the rows of the matrix of coordinates */
  for (i = 0; i < nrow; i++) {
    
    /* begin the secondary loop over the columns of the matrix of coordinates */
    for (j = 0; j < ncol; j++) {
      d[i] += pow(x1[nrow * j + i] - x2[j], 2);
    }
    
    /* matrix indexes beging at 0 in C++, while in R it starts at 1 */
    idx -= 1;
    
    /* take the squared root and replace the values in the distance matrix */
    dm(i, idx) = pow(d[i], 0.5);
  }
  return (dm);
}
/* End! */
