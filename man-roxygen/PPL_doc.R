#  Template documentation for the parameters of the family of PPL objective functions
###############################################################################################################
#' @param lags Integer value, the number of lag-distance classes. Alternatively, a vector of numeric values 
#' with the lower and upper bounds of each lag-distance class, the lowest value being larger than zero 
#' (e.g. 0.0001). Defaults to \code{lags = 7}.
#' 
#' @param lags.type Character value, the type of lag-distance classes, with options \code{"equidistant"} and
#' \code{"exponential"}. Defaults to \code{lags.type = "exponential"}.
#'
#' @param lags.base Numeric value, base of the exponential expression used to create exponentially spaced 
#' lag-distance classes. Used only when \code{lags.type = "exponential"}. Defaults to \code{lags.base = 2}.
#'
#' @param cutoff Numeric value, the maximum distance up to which lag-distance classes are created. Used only
#' when \code{lags} is an integer value. If missing, it is set to be equal to the length of the diagonal of 
#' the rectagle with sides \code{x.max} and \code{y.max} as defined in \code{\link[spsann]{scheduleSPSANN}}.
#'
#' @param criterion Character value, the feature used to describe the energy state of the system 
#' configuration, with options \code{"minimum"} and \code{"distribution"}. Defaults to 
#' \code{objective = "distribution"}.
#'
#' @param distri Numeric vector, the distribution of points or point-pairs per lag-distance class that should 
#' be attained at the end of the optimization. Used only when \code{criterion = "distribution"}. Defaults to 
#' a uniform distribution.
#' 
#' @param pairs Logical value. Should the sample configuration be optimized regarding the number of 
#' point-pairs per lag-distance class? Defaults to \code{pairs = FALSE}.
#' 
#' @concept variogram
