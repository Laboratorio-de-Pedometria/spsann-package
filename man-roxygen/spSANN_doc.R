#  Template documentation for spatial simulated annealing
################################################################################
#' @param schedule List with 10 named sub-arguments defining the control
#' parameters of the cooling schedule. See \code{scheduleSPSANN}.
#' 
#' @param plotit Logical for plotting the optimization results. This includes 
#' a) the progress of the objective function values and acceptance 
#' probabilities, and b) the original points, the perturbed points and the 
#' progress of the maximum perturbation in the x- and y-coordinates. The plots 
#' are updated at each 10 iterations. Defaults to \code{plotit = FALSE}.
#' 
#' @param boundary SpatialPolygon. The boundary of the spatial domain. 
#' If missing, it is estimated from \code{candi}.
#' 
#' @param progress Logical for printing a progress bar. Defaults to 
#' \code{progress = TRUE}.
#' 
#' @param verbose Logical for printing messages about the progress of the
#' optimization. Defaults to \code{verbose = FALSE}.
#' 
#' @param track Logical value. Should the evolution of the energy state and 
#' acceptance probability be recorded and returned with the result? If 
#' \code{track = FALSE} (the default), only the starting and ending energy state
#' values are returned with the result.
#' 
#' @keywords spatial optimize iteration
#' 
#' @concept simulated annealing

