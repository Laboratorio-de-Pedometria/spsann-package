#  Template documentation for spatial simulated annealing
################################################################################
#' @param schedule List with 10 named sub-arguments defining the control
#' parameters of the cooling schedule. See \code{scheduleSPSANN}.
#' 
#' @param plotit Logical for plotting the optimization results, including
#' a) the progress of the objective function, and b) the starting (gray) and
#' current system configuration (black), and the maximum jitter in the 
#' x- and y-coordinates. The plots are updated at each 10 jitters. Defaults to
#' \code{plotit = FALSE}.
#' 
#' @param boundary SpatialPolygon defining the boundary of the spatial domain. 
#' It is estimated from \code{candi} if missing and \code{plotit = TRUE}.
#' 
#' @param progress Logical for printing a progress bar. Defaults to 
#' \code{progress = TRUE}.
#' 
#' @param verbose Logical for printing messages about the progress of the
#' optimization. Defaults to \code{verbose = FALSE}.
#' 
#' @param track Logical value. Should the evolution of the energy state be 
#' recorded and returned with the result? If \code{track = FALSE} (the default),
#' only the starting and ending energy states are returned with the result.
#' 
#' @keywords spatial optimize iteration
#' 
#' @concept simulated annealing

