#  Template documentation for spatial simulated annealing
###############################################################################################################
#' @param schedule List with 11 named sub-arguments defining the control parameters of the cooling schedule. 
#' See \code{\link[spsann]{scheduleSPSANN}}.
#' 
#' @param plotit Logical for plotting the optimization results, including a) the progress of the objective
#' function, and b) the starting (gray) and current system configuration (black), and the maximum jitter in 
#' the x- and y-coordinates. The plots are updated at each 10 jitters. Defaults to \code{plotit = FALSE}.
#' 
#' @param boundary SpatialPolygon defining the boundary of the spatial domain. If missing and 
#' \code{plotit = TRUE}, \code{boundary} is estimated from \code{candi}.
#' 
#' @param progress Type of progress bar that should be used, with options \code{"txt"}, for a text progress 
#' bar in the R console, \code{"tk"}, to put up a Tk progress bar widget, and \code{NULL} to omit the progress
#' bar. A Tk progress bar widget is useful when using parallel processors. Defaults to 
#' \code{progress = "txt"}.
#' 
#' @param verbose Logical for printing messages about the progress of the optimization. Defaults to 
#' \code{verbose = FALSE}.
#' 
#' @param track Logical value. Should the evolution of the energy state be recorded and returned with the 
#' result? If \code{track = FALSE} (the default), only the starting and ending energy states are returned with
#' the results.
#' 
#' @keywords spatial optimize iteration
#' 
#' @concept simulated annealing

