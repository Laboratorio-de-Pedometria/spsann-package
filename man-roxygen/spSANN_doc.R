#  Template documentation for spatial simulated annealing
###############################################################################################################
#' @param points Integer value, integer vector, data frame or matrix, or list.
#' \itemize{
#' \item Integer value. The number of points. These points will be randomly sampled from \code{candi} to form
#' the starting sample configuration.
#' \item Integer vector. The row indexes of \code{candi} that correspond to the points that form the starting
#' sample configuration. The length of the vector defines the number of points.
#' \item Data frame or matrix. An object with three columns in the following order: \code{[, "id"]}, the 
#' row indexes of \code{candi} that correspond to each point, \code{[, "x"]}, the projected x-coordinates, and
#' \code{[, "y"]}, the projected y-coordinates.
#' \item List. An object with two named sub-arguments: \code{fixed}, a data frame or matrix with the projected
#' x- and y-coordinates of the existing sample configuration -- kept fixed during the optimization --, and 
#' \code{free}, an integer value defining the number of points that should be added to the existing sample
#' configuration -- free to move during the optimization.
#' }
#' 
#' @param candi Data frame or matrix with the candidate locations for the jittered points. \code{candi} must 
#' have two columns in the following order: \code{[, "x"]}, the projected x-coordinates, and \code{[, "y"]}, 
#' the projected y-coordinates.
#' 
#' @param schedule List with 11 named sub-arguments defining the control parameters of the cooling schedule. 
#' See \code{\link[spsann]{scheduleSPSANN}}.
#' 
#' @param plotit (Optional) Logical for plotting the optimization results, including a) the progress of the
#' objective function, and b) the starting (gray circles) and current sample configuration (black dots), and 
#' the maximum jitter in the x- and y-coordinates. The plots are updated at each 10 jitters. When adding 
#' points to an existing sample configuration, fixed points are indicated using black crosses. Defaults to 
#' \code{plotit = FALSE}.
#' 
#' @param boundary (Optional) SpatialPolygon defining the boundary of the spatial domain. If missing and 
#' \code{plotit = TRUE}, \code{boundary} is estimated from \code{candi}.
#' 
#' @param progress (Optional) Type of progress bar that should be used, with options \code{"txt"}, for a text
#' progress bar in the R console, \code{"tk"}, to put up a Tk progress bar widget, and \code{NULL} to omit the
#' progress bar. A Tk progress bar widget is useful when using parallel processors. Defaults to 
#' \code{progress = "txt"}.
#' 
#' @param verbose (Optional) Logical for printing messages about the progress of the optimization. Defaults to 
#' \code{verbose = FALSE}.
#' 
#' @param track (Optional) Logical value. Should the evolution of the energy state be recorded and returned 
#' along with the result? If \code{track = FALSE} (the default), only the starting and ending energy states 
#' are returned along with the results.
#' 
#' @keywords spatial optimize iteration
#' 
#' @concept simulated annealing

