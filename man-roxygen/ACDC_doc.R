#  Template documentation for the family of ACDC objective functions
################################################################################
#' @param strata.type Character value setting the type of stratification that 
#' should be used to create sampling strata (or factor levels) for the numeric 
#' covariates. Available options are \code{"area"} for equal-area, and
#' \code{"range"} for equal-range. Defaults to \code{strata.type = "area"}. See
#' \sQuote{Details} for more information.
#' 
#' @param covars Data frame or matrix with the covariates in the columns.
#'
#' @param use.coords Logical value. Should the geographic coordinates be used as
#' covariates? Defaults to \code{use.coords = FALSE}.
#'
#' @references
#' Cramér, H. \emph{Mathematical methods of statistics}. Princeton: Princeton 
#' University Press, p. 575, 1946.
#'
#' Everitt, B. S. \emph{The Cambridge dictionary of statistics}. Cambridge: 
#' Cambridge University Press, p. 432, 2006.
#'
#' Hyndman, R. J.; Fan, Y. Sample quantiles in statistical packages. \emph{The 
#' American Statistician}, v. 50, p. 361-365, 1996.
#'
#' Minasny, B.; McBratney, A. B. A conditioned Latin hypercube method for
#' sampling in the presence of ancillary information. \emph{Computers &
#' Geosciences}, v. 32, p. 1378-1388, 2006.
#'
#' Minasny, B.; McBratney, A. B. Conditioned Latin Hypercube Sampling for
#' calibrating soil sensor data to soil properties. Chapter 9. Viscarra Rossel,
#' R. A.; McBratney, A. B.; Minasny, B. (Eds.) \emph{Proximal Soil Sensing}.
#' Amsterdam: Springer, p. 111-119, 2010.
#'
#' Mulder, V. L.; de Bruin, S.; Schaepman, M. E. Representing major soil
#' variability at regional scale by constrained Latin hypercube sampling of
#' remote sensing data. \emph{International Journal of Applied Earth Observation
#' and Geoinformation}, v. 21, p. 301-310, 2013.
#'
#' Roudier, P.; Beaudette, D.; Hewitt, A. A conditioned Latin hypercube sampling
#' algorithm incorporating operational constraints. \emph{5th Global Workshop on
#' Digital Soil Mapping}. Sydney, p. 227-231, 2012.
#'
#' @note
#' This function was derive with modifications from the method known as the 
#' \emph{conditioned Latin Hypercube sampling} originally proposed by Minasny 
#' and McBratney (2006).


