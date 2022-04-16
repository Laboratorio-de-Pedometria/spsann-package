#  Template documentation for parameters of the family of ACDC objective functions
####################################################################################################
#' @param strata.type (Optional) Character value setting the type of stratification that should be
#' used to create the marginal sampling strata (or factor levels) for the numerical covariates. Two
#' options are available:
#'  * `"area"` (Default) Equal-area marginal sampling strata.
#'  * `"range"` Equal-range marginal sampling strata.
#'
#' The first option (`"area"`) is equivalent to drawing the frequency histogram of the numerical
#' covariates with bins of variable width but equal area. The second, however, would results in a
#' frequency histogram with bins of equal width but variable area such as when using
#' [graphics::hist()] with the default options. Marginal sampling strata of equal-area will include
#' virtually the same number of individual covariate cells, while equal-range strata aim for the
#' same number of individual covariate values.
#'
#' @param covars Data frame or matrix with the spatially exhaustive covariates in the columns.
#'
#' @param use.coords (Optional) Logical value. Should the projected spatial x- and y-coordinates
#' be used as spatially exhaustive covariates? Defaults to `use.coords = FALSE`.
#'
#' @concept spatial trend
