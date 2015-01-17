#  Template documentation for spatial jittering
################################################################################
#' @param points Integer value, integer vector, data frame or matrix. If 
#' \code{points} is an integer value, it defines the number of points that 
#' should be sampled from \code{candi}. If \code{points} is an integer vector, 
#' it contains the row indexes of \code{candi} that correspond to the points. 
#' If \code{points} is a data frame or matrix, it must have three
#' columns: \code{[, 1]} the identification of each point (1, 2, ..., n), 
#' \code{[, 2]} the x coordinates, and \code{[, 3]} the y coordinates. The
#' coordinates must be projected. If a data frame or matrix is used, 
#' \code{points} must be a subset of \code{candi}.
#'
#' @param candi Data frame or matrix with the candidate locations for the
#' sample points. The coordinates must be projected.
#'
#' @param x.max,x.min,y.max,y.min Numeric value. The minimum and maximum 
#' quantity of random noise to be added to the x and y coordinates. The minimum
#' quantity must be equal to, at least, the minimum distance between two
#' neighboring candidate locations. The units are the same as of the coordinates.
