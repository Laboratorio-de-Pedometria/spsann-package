#  Template documentation for parameters of the family of ACDC objective functions
###############################################################################################################
#' @param strata.type (Optional) Character value setting the type of stratification that should be used to 
#' create the marginal sampling strata (or factor levels) for the numeric covariates. Available options are
#' \code{"area"}, for equal-area, and \code{"range"}, for equal-range. Defaults to \code{strata.type = "area"}.
#' 
#' @param covars Data frame or matrix with the covariates in the columns.
#'
#' @param use.coords (Optional) Logical value. Should the spatial x- and y-coordinates be used as covariates? 
#' Defaults to \code{use.coords = FALSE}.
#'
#' @concept spatial trend
