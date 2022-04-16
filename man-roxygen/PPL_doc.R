#  Template documentation for the parameters of the family of PPL objective functions
####################################################################################################
#' @param lags Integer value, the number of lag-distance classes. Alternatively, a vector of numeric
#' values with the lower and upper bounds of each lag-distance class, the lowest value being larger
#' than zero (e.g. 0.0001). Defaults to `lags = 7`.
#'
#' @param lags.type Character value, the type of lag-distance classes, with options `"equidistant"`
#' and `"exponential"`. Defaults to `lags.type = "exponential"`.
#'
#' @param lags.base Numeric value, base of the exponential expression used to create exponentially
#' spaced lag-distance classes. Used only when `lags.type = "exponential"`. Defaults to
#' `lags.base = 2`.
#'
#' @param cutoff Numeric value, the maximum distance up to which lag-distance classes are created.
#' Used only when `lags` is an integer value. If missing, it is set to be equal to the length of
#' the diagonal of the rectangle with sides `x.max` and `y.max` as defined in
#' [spsann::scheduleSPSANN()].
#'
#' @param criterion Character value, the feature used to describe the energy state of the system
#' configuration, with options `"minimum"` and `"distribution"`. Defaults to
#' `objective = "distribution"`.
#'
#' @param distri Numeric vector, the distribution of points or point-pairs per lag-distance class
#' that should be attained at the end of the optimization. Used only when
#' `criterion = "distribution"`. Defaults to a uniform distribution.
#'
#' @param pairs Logical value. Should the sample configuration be optimized regarding the number of
#' point-pairs per lag-distance class? Defaults to `pairs = FALSE`.
#'
#' @concept variogram
