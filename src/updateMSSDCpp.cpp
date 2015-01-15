/*******************************************************************************
Update the matrix of distances
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
  NumericVector d(nrow, 0.0000);
  
  for (i = 0; i < nrow; i++) { /* loop over the rows */    
    for (j = 0; j < ncol; j++) { /* loop over the columns */
      d[i] += pow(x1[nrow * j + i] - x2[j], 2);
    }
    /* take the squared root and replace the values in the distance matrix */
    dm(i, idx - 1) = pow(d[i], 0.5);
  }
  return (dm);
}
// # Testing
// require(SpatialTools)
// x1 <- matrix(rep(1, 6), ncol = 2); x1
// x1[2, ] <- c(2, 2); x1
// b <- x1[c(1, 3), ]; b
// dm <- SpatialTools::dist2(x1, b); dm
// idx <- as.integer(2); idx
// x2 <- matrix(x1[idx, ], nrow = 1); x2
// .updateMSSDCpp(x1, x2, dm, idx)
