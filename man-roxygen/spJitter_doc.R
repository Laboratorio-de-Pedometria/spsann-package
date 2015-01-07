#  Template documentation for spatial jittering
################################################################################
#' @param points Integer value setting the number of points. Alternatively, a
#' data.frame or matrix with three columns: 1) the identification of each point,
#' 2) the x coordinates, and 3) the y coordinates. The coordinates must be
#' projected. If a data.frame or matrix is used, \code{points} must be a subset
#' of \code{candidates}. See \sQuote{Details} for more information.
#' 
#' @param candidates data.frame or matrix with the candidate locations for the 
#' sample points. See \sQuote{Details} for more information.
#' 
#' @param x.max,x.min,y.max,y.min The minimum and maximum quantity of random 
#' noise to be added to the x and y coordinates. The minimum quantity should be
#' equal to, at least, the minimum distance between two neighboring candidate 
#' locations. The units are the same as of the coordinates. See \sQuote{Details}
#' for more information.
