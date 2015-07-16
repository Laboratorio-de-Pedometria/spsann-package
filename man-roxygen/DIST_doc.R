#  Template documentation for the family of DIST objective functions
################################################################################
#' @section Marginal distribution of covariates:
#' Reproducing the marginal distribution of the numeric covariates depends upon
#' the definition of marginal sampling strata. These marginal sampling strata 
#' are also used to define the factor levels of all numeric covariates that  
#' are passed together with factor covariates.
#' 
#' Two types of marginal sampling strata can be used. \emph{Equal-area} 
#' sampling strata are defined using the sample quantiles estimated with 
#' \code{\link[stats]{quantile}} using a discontinuous function 
#' (\code{type = 3}). This is to avoid creating breakpoints that do not occur 
#' in the population of existing covariate values.
#' 
#' The function \code{\link[stats]{quantile}} commonly produces repeated 
#' break points. A break point will always be repeated if that value has a 
#' relatively high frequency in the population of covariate values. The number 
#' of repeated break points increases with the number of marginal sampling 
#' strata. Only unique break points are used to create marginal sampling strata.
#' 
#' \emph{Equal-range} sampling strata are defined by breaking the range of 
#' covariate values into pieces of equal size. This method usually creates
#' break points that do not occur in the population of existing covariate 
#' values. Such break points are replaced by the nearest existing covariate 
#' value identified using Euclidean distances.
#' 
#' Both stratification methods can produce marginal sampling strata that cover 
#' a range of values that do not exist in the population of covariate values. 
#' Any empty marginal sampling strata is merged with the closest non-empty 
#' marginal sampling strata. These are identified using Euclidean distances.
#' 
#' The approaches used to define the marginal sampling strata result in each 
#' numeric covariate having a different number of marginal sampling strata, 
#' some of them with different area/size. Because the goal is to have a sample 
#' that reproduces the marginal distribution of the covariate, each marginal 
#' sampling strata will have a different number of sample points. The wanted 
#' distribution of the number of sample points per marginal strata is estimated 
#' empirically as the proportion of points in the population of existing 
#' covariate values that fall in each marginal sampling strata. 


