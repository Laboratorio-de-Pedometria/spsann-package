#  Template documentation for spatial simulated annealing
####################################################################################################
#' @param points Integer value, integer vector, data frame or matrix, or list. The number of
#' sampling points (sample size) or the starting sample configuration. Four options are available:
#'  * Integer value. The required number of sampling points (sample size). The sample configuration
#'    used to start the optimization will consist of grid cell centres of `candi` selected using
#'    simple random sampling, i.e. [base::sample()] with `x = 1:nrow(candi)` and `size = points`.
#'  * Integer vector. A set of row indexes between one (1) and `nrow(candi)` identifying the grid
#'    cell centres of `candi` that should be used to form the starting sample configuration for the
#'    optimization. The length of the integer vector, `length(points)`, is the sample size.
#'  * Data frame or matrix. The Cartesian x- and y-coordinates of the starting sample configuration.
#'  * List. An object with two named sub-arguments:
#'    * `fixed` An integer vector, data frame or matrix specifying an existing sample configuration
#'      (see options above). This sample configuration is kept as-is (fixed) during the entire
#'      optimization, being used only to compute the objetive function values.
#'    * `free` An integer value, integer vector, data frame or matrix (see options above) specifying
#'      the (number of) sampling points that should be added to the existing sample configuration.
#'      These new sampling points are free to be moved around (jittered) during the optimization.
#'
#'  Most users will want to set an integer value simply specifying the required sample size. Using
#'  an integer vector or data frame (or matrix) will generally be useful to users willing to
#'  evaluate starting sample configurations, test strategies to speed up the optimization, and
#'  fine-tune or thin an existing sample configuration. Finally, a list, will generally be used to
#'  augment a possibly already existing, real-world sample configuration or fine-tune only a subset
#'  of the existing sampling points.
#'
#' @param candi Data frame or matrix with the candidate locations for the jittered samples. `candi`
#' must have two columns in the following order: `[, "x"]`, the Cartesian x-coordinates, and
#' `[, "y"]`, the Cartesian y-coordinates.
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
