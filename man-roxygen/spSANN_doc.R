#  Template documentation for spatial simulated annealing
####################################################################################################
#' @param points Integer value, integer vector, data frame or matrix, or list.
#'  * Integer value. The required number of sampling points (sample size). The sample configuration
#'    used to start the optimization will consist of grid cell centres of `candi` selected using
#'    simple random sampling, i.e. [base::sample()] with `x = 1:nrow(candi)` and `size = points`.
#'  * Integer vector. A set of row indexes between one (1) and `nrow(candi)` identifying the grid
#'    cell centres of `candi` that should be used to form the starting sample configuration for the
#'    optimization. The length of the integer vector, `length(points)`, defines the number of
#'    sampling points (sample size).
#' 
#'  * Data frame or matrix. An object with three columns in the following order: `[, "id"]`, the
#'    row indexes of `candi` that correspond to each sample, `[, "x"]`, the projected x-coordinates,
#'    and `[, "y"]`, the projected y-coordinates.
#'  * List. An object with two named arguments:
#'    * `fixed`: a data frame or matrix with the projected x- and y-coordinates of the existing
#'      sample configuration -- kept fixed during the optimization --, and
#'    * `free`: an integer value defining the number of samples that should be added to the existing
#'      sample configuration -- free to move during the optimization.
#'
#' @param candi Data frame or matrix with the candidate locations for the jittered samples. `candi`
#' must have two columns in the following order: `[, "x"]`, the projected x-coordinates, and
#' `[, "y"]`, the projected y-coordinates.
#'
#' @param schedule List with named sub-arguments defining the control parameters of the cooling
#' schedule. See [spsann::scheduleSPSANN()].
#'
#' @param plotit (Optional) Logical for plotting the optimization results, including a) the progress
#' of the objective function, and b) the starting (gray circles) and current sample configuration
#' (black dots), and the maximum jitter in the x- and y-coordinates. The plots are updated at each
#' 10 jitters. When adding samples to an existing sample configuration, fixed samples are indicated
#' using black crosses. Defaults to `plotit = FALSE`.
#'
#' @param boundary (Optional) SpatialPolygon defining the boundary of the spatial domain. If
#' missing and `plotit = TRUE`, `boundary` is estimated from `candi`.
#'
#' @param progress (Optional) Type of progress bar that should be used, with options `"txt"`, for a
#' text progress bar in the R console, `"tk"`, to put up a Tk progress bar widget, and `NULL` to
#' omit the progress bar. A Tk progress bar widget is useful when using parallel processors.
#' Defaults to `progress = "txt"`.
#'
#' @param verbose (Optional) Logical for printing messages about the progress of the optimization.
#' Defaults to `verbose = FALSE`.
#'
#' @param track (Optional) Logical value. Should the evolution of the energy state be recorded and
#' returned along with the result? If `track = FALSE` (the default), only the starting and ending
#' energy states are returned along with the results.
#'
#' @keywords spatial optimize iteration
#'
#' @concept simulated annealing
