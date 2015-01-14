#  Template documentation for optimPPL()
################################################################################
#' @param lags Integer. The number of lag distance classes. Alternatively, a
#' vector of numeric values with the lower and upper limits of each lag
#' distance class. The lowest value must be larger than zero, e.g. 0.0001.
#' Defaults to \code{lags = 7}.
#'
#' @param lags.type Character. The type of lag distance classes. Available
#' options are \code{"equidistant"} and \code{"exponential"}. Defaults to
#' \code{lags.type = "exponential"}. See \sQuote{Details} for more information.
#'
#' @param lags.base Numeric. Base of the exponential expression used to
#' create the exponential lag distance classes. Defaults to
#' \code{lags.base = 2}. See \sQuote{Details} for more information.
#'
#' @param cutoff Numeric. The maximum distance up to which lag distance classes
#' are created. Used only when \code{lags} is an integer. See \sQuote{Details}
#' for more information.
#'
#' @param criterion Character. The measure to be used to describe the energy
#' state of the current system configuration. Available options are
#' \code{"minimum"} and \code{"distribution"}. Defaults to
#' \code{objective = "minimum"}. See \sQuote{Details} for more information.
#'
#' @param pre.distri Numeric vector. The pre-specified distribution of points
#' or point-pair with which the observed counts of points or point-pairs per
#' lag distance class is compared. Used only when
#' \code{criterion = "distribution"}. Defaults to a uniform distribution. See
#' \sQuote{Details} for more information.
