#  Template documentation for multi-objective optimization problems (MOOP)
################################################################################
#' @param nadir List with named sub-arguments. Three options are available: 
#' 1) \code{sim} -- the number of simulations that should be used to estimate 
#' the nadir point, and \code{seeds} -- vector defining the random seeds for
#' each simulation; 2) \code{user} -- a list of user-defined nadir values named 
#' after the respective objective function to which they apply; 3) \code{abs} 
#' -- logical for calculating the nadir point internally (experimental).
#'
#' @param weights List with named sub-arguments. The weights assigned to each 
#' one of the objective functions that form the multi-objective optimization
#' problem (MOOP). They must be named after the respective objective function 
#' to which they apply. The weights must be equal to or larger than 0 and sum 
#' to 1. The default option gives equal weights to all objective functions.
#'
#' @param utopia List with named sub-arguments. Two options are available: 1) 
#' \code{user} -- a list of user-defined values named after the respective 
#' objective function to which they apply; 2) \code{abs} -- logical for 
#' calculating the utopia point internally (experimental).
#' 
#' @section Multi-objective optimization:
#' A method of solving a multi-objective optimization problem is to aggregate 
#' the objective functions into a single \emph{utility function}. In the
#' \pkg{spsann} package, the aggregation is performed using the \emph{weighted 
#' sum method}, which incorporates in the weights the preferences of the user 
#' regarding the relative importance of each objective function.
#' 
#' The weighted sum method is affected by the relative magnitude of the 
#' different function values. The objective functions implemented in the
#' \pkg{spsann} package have different units and orders of magnitude. The 
#' consequence is that the objective function with the largest values will have 
#' a numerical dominance in the optimization. In other words, the weights will 
#' not express the true preferences of the user, and the meaning of the utility 
#' function becomes unclear.
#' 
#' A solution to avoid the numerical dominance is to transform the objective
#' functions so that they are constrained to the same approximate range of 
#' values. Several function-transformation methods can be used and the 
#' \pkg{spsann} offers a few of them. The \emph{upper-lower-bound approach}
#' requires the user to inform the maximum (nadir point) and minimum (utopia
#' point) absolute function values. The resulting function values will always 
#' range between 0 and 1.
#' 
#' Using the \emph{upper-bound approach} requires the user to inform only the
#' nadir point, while the utopia point is set to zero. The upper-bound approach
#' for transformation aims at equalizing only the upper bounds of the objective 
#' functions. The resulting function values will always be smaller than or equal
#' to 1.
#' 
#' Sometimes, the absolute maximum and minimum values of an objective function 
#' can be calculated exactly. This seems not to be the case of the objective 
#' functions implemented in the \pkg{spsann} package. If the user is 
#' uncomfortable with informing the nadir and utopia points, there is the option
#' for using \emph{numerical simulations}. It consists in computing the 
#' function value for many random sample configurations. The mean function 
#' value is used to set the nadir point, while the the utopia point is set to
#' zero. This approach is similar to the upper-bound approach, but the function
#' values will have the same orders of magnitude only at the starting point of 
#' the optimization. Function values larger than one are likely to occur during 
#' the optimization. We recommend the user to avoid this approach whenever 
#' possible because the effect of the starting point on the optimization as a 
#' whole usually is insignificant or arbitrary.
#' 
#' The \emph{upper-lower-bound approach} with the \emph{Pareto maximum and 
#' minimum values} is the most elegant solution to transform the objective 
#' functions. However, it is the most time consuming. It works as follows:
#' 
#' \enumerate{
#'   \item Optimize a sample configuration with respect to each objective
#'   function that composes the MOOP;
#'   \item Compute the function value of every objective function for every
#'   optimized sample configuration;
#'   \item Record the maximum and minimum absolute function values computed for 
#'   each objective function--these are the Pareto maximum and minimum.
#' }
#' 
#' For example, consider that a MOOP is composed of two objective functions (A 
#' and B). The minimum absolute value for function A is obtained when the sample
#' configuration is optimized with respect to function A. This is the Pareto
#' minimum of function A. Consequently, the maximum absolute value for function
#' A is obtained when the sample configuration is optimized regarding function
#' B. This is the Pareto maximum of function A. The same logic applies for 
#' function B.
#'
#' @references
#' Arora, J. \emph{Introduction to optimum design}. Waltham: Academic Press, p. 
#' 896, 2011.
#'
#' Marler, R. T.; Arora, J. S. Survey of multi-objective optimization methods 
#' for engineering. \emph{Structural and Multidisciplinary Optimization}, v. 26,
#' p. 369-395, 2004.
#' 
#' Marler, R. T.; Arora, J. S. Function-transformation methods for 
#' multi-objective optimization. \emph{Engineering Optimization}, v. 37, p. 
#' 551-570, 2005.
#'
#' Marler, R. T.; Arora, J. S. The weighted sum method for multi-objective 
#' optimization: new insights. \emph{Structural and Multidisciplinary 
#' Optimization}, v. 41, p. 853-862, 2009.
