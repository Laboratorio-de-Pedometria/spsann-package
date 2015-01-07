/*******************************************************************************
Calculate the mean squared shortest distance
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
