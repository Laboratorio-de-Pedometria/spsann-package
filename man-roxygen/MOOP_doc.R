#  Template documentation for multi-objective optimization problems (MOOP)
################################################################################
#'
#' @param nadir List with named sub-arguments. Three options are available: 
#' 1) \code{sim} -- the number of random realizations to estimate the nadir 
#' point, and \code{seeds} -- vector defining the random seeds for each of the
#' realizations; 2) \code{user} -- a list of user-defined values named after 
#' the objective function to which they apply; 3) \code{abs} -- logical for 
#' calculating the nadir point internally.
#'
#' @param weights List with named sub-arguments. The weights assigned to each of
#' the objective functions combined to form the multi-objective optimization
#' problem (MOOP). They must be named after the objective function to which
#' they apply. The weights must be larger than 0 and sum to 1. The default
#' option gives equal weights to all objective functions.
#'
#' @param utopia List with two named sub-arguments: \code{user} -- a list of
#' user-defined values named after the objective function to which they apply, 
#' and \code{abs} -- logical for calculating the utopia point internally. 
#' Defaults to \code{user = NULL} and \code{abs = NULL}.
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


