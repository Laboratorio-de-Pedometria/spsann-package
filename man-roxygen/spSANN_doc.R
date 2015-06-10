#  Template documentation for spatial simulated annealing
################################################################################
#' @param iterations Integer. The maximum number of iterations that should be
#' used for the optimization.
#' 
#' @param acceptance List with two named sub-arguments: \code{initial} -- 
#' numeric value between 0 and 1 defining the initial acceptance probability, 
#' and \code{cooling} -- a numeric value defining the exponential factor by
#' witch the acceptance probability decreases at each iteration. Defaults to 
#' \code{acceptance = list(initial = 0.99, cooling = iterations / 10)}.
#' 
#' @param stopping List with one named sub-argument: \code{max.count} -- 
#' integer value defining the maximum allowable number of iterations without 
#' improvement of the objective function value. Defaults to 
#' \code{stopping = list(max.count = iterations / 10)}. More options may be
#' included in the future.
#' 
#' @param plotit Logical for plotting the optimization results. This includes 
#' a) the progress of the objective function values and acceptance 
#' probabilities, and b) the original points, the perturbed points and the 
#' progress of the maximum perturbation in the x- and y-coordinates. The plots 
#' are updated at each 10 iterations. Defaults to \code{plotit = TRUE}.
#' 
#' @param boundary SpatialPolygon. The boundary of the spatial domain. 
#' Mandatory if \code{plotit = TRUE}.
#' 
#' @param progress Logical for printing a progress bar. Defaults to 
#' \code{progress = TRUE}.
#' 
#' @param verbose Logical for printing messages about the progress of the
#' optimization.
#' 
#' @param greedy Logical value. Should the optimization be done using a greedy
#' algorithm, that is, accepting only better system configurations? Defaults
#' to \code{greedy = FALSE}.
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
#' @keywords spatial optimize
#' 
#' @concept simulated annealing
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
