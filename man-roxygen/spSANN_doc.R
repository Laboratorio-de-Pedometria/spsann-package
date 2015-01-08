#  Template documentation for spatial simulated annealing
################################################################################
#' @param iterations Integer. The maximum number of iterations that should be
#' used for the optimization.
#' 
#' @param acceptance List with two sub-arguments: \code{initial} and 
#' \code{cooling}. \code{initial} is a numeric value between 0 and 1 defining
#' the initial acceptance probability. Defaults to \code{initial = 0.99}.
#' \code{cooling} is a numeric value defining the exponential factor by with 
#' the acceptance probability decreases at each iteration. Defaults to 
#' \code{cooling = iterations / 10}.
#' 
#' @param stopping List with one sub-argument: \code{max.count}. 
#' \code{max.count} is an integer value defining the maximum allowable number 
#' of iterations without improvement of the objective function value. Defaults 
#' to \code{max.count = iterations / 10}.
#' 
#' @param plotit Logical for ploting the optimization results. This includes 
#' a) the progress of the objective function values and acceptance 
#' probabilities, and b) the original points, the perturbed points and the 
#' progress of the maximum perturbation in the x and y coordinates. The plots 
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
#' @keywords spatial optimize
