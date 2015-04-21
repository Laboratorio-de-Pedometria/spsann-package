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
#' @keywords spatial optimize
#' 
#' @concept simulated annealing
