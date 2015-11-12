#  Template documentation for the family of DIST objective functions
################################################################################
#' @section Marginal distribution of covariates:
#' Reproducing the marginal distribution of the numeric covariates depends upon
#' the definition of marginal sampling strata. These marginal sampling strata 
#' are also used to define the factor levels of all numeric covariates that  
#' are passed together with factor covariates. Two types of marginal sampling 
#' strata can be used: \emph{equal-area} and \emph{equal-range}.
#' 
#' \subsection{Equal-area}{
#' \emph{Equal-area} marginal sampling strata are defined using the sample 
#' quantiles estimated with \code{\link[stats]{quantile}} using a discontinuous 
#' function (\code{type = 3}). This is to avoid creating breakpoints that do 
#' not occur in the population of existing covariate values.
#' 
#' Depending on the level of discretization of the covariate values, 
#' \code{\link[stats]{quantile}} produces repeated breakpoints. A breakpoint 
#' will be repeated if that value has a relatively high frequency in the 
#' population of covariate values. The number of repeated breakpoints increases 
#' with the number of marginal sampling strata. Repeated breakpoints result in
#' empty marginal sampling strata. To avoid this, only the unique breakpoints 
#' are used.
#' }
#' \subsection{Equal-range}{
#' \emph{Equal-range} marginal sampling strata are defined by breaking the range
#' of covariate values into pieces of equal size. Depending on the level of 
#' discretization of the covariate values, this method creates breakpoints that
#' do not occur in the population of existing covariate values. Such breakpoints
#' are replaced by the nearest existing covariate value identified using 
#' Euclidean distances.
#' 
#' Like the equal-area method, the equal-range method can produce empty marginal
#' sampling strata. The solution used here is to merge any empty marginal 
#' sampling strata with the closest non-empty marginal sampling strata. This is
#' identified using Euclidean distances as well.
#' }
#' 
#' The approaches used to define the marginal sampling strata result in each 
#' numeric covariate having a different number of marginal sampling strata, 
#' some of them with different area/size. Because the goal is to have a sample 
#' that reproduces the marginal distribution of the covariate, each marginal 
#' sampling strata will have a different number of sample points. The wanted 
#' distribution of the number of sample points per marginal strata is estimated 
#' empirically as the proportion of points in the population of existing 
#' covariate values that fall in each marginal sampling strata.
