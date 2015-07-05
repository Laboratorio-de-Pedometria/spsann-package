#  Template documentation for spatial jittering
################################################################################
#' @param points Integer value, integer vector, data frame or matrix. If 
#' \code{points} is an integer value, it defines the number of points that 
#' should be randomly sampled from \code{candi} to form the starting system
#' configuration. If \code{points} is a vector of integer values, it contains 
#' the row indexes of \code{candi} that correspond to the points that form the
#' starting system configuration. If \code{points} is a data frame or matrix, 
#' it must have three columns in the following order: \code{[, "id"]} the 
#' row indexes of \code{candi} that correspond to each point, \code{[, "x"]} 
#' the projected x-coordinates, and \code{[, "y"]} the projected y-coordinates.
#' Note that in the later case, \code{points} must be a subset of \code{candi}.
#'
#' @param candi Data frame or matrix with the candidate locations for the
#' perturbed points. \code{candi} must have two columns in the following 
#' order: \code{[, "x"]} the projected x-coordinates, and \code{[, "y"]} the 
#' projected y-coordinates.
#'
#' @param x.max,x.min,y.max,y.min Numeric value. The minimum and maximum 
#' quantity of random noise to be added to the projected x- and y-coordinates.
#' The minimum quantity should be equal to, at least, the minimum distance
#' between two neighbouring candidate locations. The units are the same as of 
#' the projected x- and y-coordinates. If missing, they are estimated from
#' \code{candi}.
