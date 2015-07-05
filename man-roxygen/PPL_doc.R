#  Template documentation for the family of PPL objective functions
################################################################################
#' @param lags Integer value. The number of lag-distance classes. Alternatively,
#' a vector of numeric values with the lower and upper limits of each 
#' lag-distance class. The lowest value must be larger than zero. Defaults to
#' \code{lags = 7}.
#' 
#' @param lags.type Character value. The type of lag-distance classes, with
#' options \code{"equidistant"} and \code{"exponential"}. Defaults to
#' \code{lags.type = "exponential"}.
#'
#' @param lags.base Numeric value. Base of the exponential expression used to
#' create exponentially spaced lag-distance classes. Used only when 
#' \code{lags.type = "exponential"}. Defaults to \code{lags.base = 2}.
#'
#' @param cutoff Numeric value. The maximum distance up to which lag-distance
#' classes are created. Used only when \code{lags} is an integer value. 
#'
#' @param criterion Character value. The feature used to describe the
#' energy state of the system configuration, with options \code{"minimum"} and
#' \code{"distribution"}. Defaults to \code{objective = "distribution"}.
#'
#' @param distri Numeric vector. The distribution of points or point-pairs per
#' lag-distance class that should be attained at the end of the optimization. 
#' Used only when \code{criterion = "distribution"}. Defaults to a uniform
#' distribution.
#' 
#' @param pairs Logical value. Should the sample configuration be optimized
#' regarding the number of point-pairs per lag-distance class? Defaults to 
#' \code{pairs = FALSE}.
#'
#' @references
#' Bresler, E.; Green, R. E. \emph{Soil parameters and sampling scheme for
#' characterizing soil hydraulic properties of a watershed}. Honolulu:
#' University of Hawaii at Manoa, p. 42, 1982.
#'
#' Marler, R. T.; Arora, J. S. Function-transformation methods for
#' multi-objective optimization. \emph{Engineering Optimization}. v. 37, p.
#' 551-570, 2005.
#' 
#' Pettitt, A. N. & McBratney, A. B. Sampling designs for estimating spatial
#' variance components. \emph{Applied Statistics}. v. 42, p. 185, 1993.
#'
#' Russo, D. Design of an optimal sampling network for estimating the variogram.
#' \emph{Soil Science Society of America Journal}. v. 48, p. 708-716, 1984.
#'
#' Truong, P. N.; Heuvelink, G. B. M.; Gosling, J. P. Web-based tool for expert
#' elicitation of the variogram. \emph{Computers and Geosciences}. v. 51, p.
#' 390-399, 2013.
#'
#' Warrick, A. W.; Myers, D. E. Optimization of sampling locations for variogram
#' calculations. \emph{Water Resources Research}. v. 23, p. 496-500, 1987.
