#  Template documentation for spatial jittering (include in functions other than spJitter)
###############################################################################################################
#' @details
#' Details about the mechanism used to generate a new sample configuration out of the current sample 
#' configuration by randomly perturbing the coordinates of a sample point are available in the help page of
#' \code{\link[spsann]{spJitter}}.
#' 
#' @note 
#' The distance between two points is computed as the Euclidean distance between them. This computation 
#' assumes that the optimization is operating in the two-dimensional Euclidean space, i.e. the coordinates of
#' the sample points and candidate locations should not be provided as latitude/longitude. \pkg{spsann} has no 
#' mechanism to check if the coordinates are projected: the user is responsible for making sure that this
#' requirement is attained.
