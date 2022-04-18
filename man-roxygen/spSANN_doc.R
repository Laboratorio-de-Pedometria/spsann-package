#  Template documentation for spatial simulated annealing
####################################################################################################
#' @param points Integer value, integer vector, data frame (or matrix), or list. The number of
#' sampling points (sample size) or the starting sample configuration. Four options are available:
#'  * Integer value. The required number of sampling points (sample size). The sample configuration
#'    used to start the optimization will consist of grid cell centres of `candi` selected using
#'    simple random sampling, i.e. [base::sample()] with `x = 1:nrow(candi)` and `size = points`.
#'  * Integer vector. A set of row indexes between one (1) and `nrow(candi)`. These row indexes
#'    identify the grid cell centres of `candi` that will form the starting sample configuration for
#'    the optimization. The length of the integer vector, `length(points)`, is the sample size.
#'  * Data frame (or matrix). The Cartesian x- and y-coordinates (in this order) of the starting
#'    sample configuration.
#'  * List. An object with two named sub-arguments:
#'    * `fixed` An integer vector or data frame (or matrix) specifying an existing sample
#'      configuration (see options above). This sample configuration is kept as-is (fixed) during
#'      the optimization and is used only to compute the objective function values.
#'    * `free` An integer value, integer vector, data frame or matrix (see options above) specifying
#'      the (number of) sampling points to add to the existing sample configuration. These new
#'      sampling points are free to be moved around (jittered) during the optimization.
#'
#'  Most users will want to set an integer value simply specifying the required sample size. Using
#'  an integer vector or data frame (or matrix) will generally be helpful to users willing to
#'  evaluate starting sample configurations, test strategies to speed up the optimization, and
#'  fine-tune or thin an existing sample configuration. Users interested in augmenting a possibly
#'  existing real-world sample configuration or fine-tuning only a subset of the existing sampling
#'  points will want to use a list.
#'
#' @param candi Data frame (or matrix). The Cartesian x- and y-coordinates (in this order) of the
#' cell centres of a spatially exhaustive, rectangular grid covering the entire spatial sampling
#' domain. The spatial sampling domain can be contiguous or composed of disjoint areas and contain
#' holes and islands. `candi` provides the set of (finite) candidate locations inside the spatial
#' sampling domain for a point jittered during the optimization. Usually, `candi` will match the
#' geometry of the spatial grid containing the prediction locations, e.g. `newdata`
#' in [gstat::krige()], `object` in [raster::predict()], and `locations` in [geoR::krige.conv()].
#'
#' @param schedule List with named sub-arguments setting the control parameters of the annealing
#' schedule. See [spsann::scheduleSPSANN()].
#'
#' @param plotit (Optional) Logical for plotting the evolution of the optimization. Plot updates
#' occur at each ten (10) spatial jitters. Defaults to `plotit = FALSE`. The plot includes two
#' panels:
#'  1. The first panel depicts the changes in the objective function value (y-axis) with the
#'    annealing schedule (x-axis). The objective function values should be high and variable at the
#'    beginning of the optimization (panel's top left). As the optimization proceeds, the objective
#'    function values should gradually transition to a monotone decreasing behaviour till they
#'    become virtually constant. The objective function values constancy suggests the end of the
#'    optimization (panel's bottom right).
#'  2. The second panel shows the starting (grey circles) and current spatial sample configuration
#'    (black dots). Black crosses indicate the fixed (existing) sampling points when a spatial
#'    sample configuration is augmented. The plot shows the starting sample configuration to assess
#'    the effects on the optimized spatial sample configuration: the latter generally should be
#'    independent of the first. The second panel also shows the maximum possible spatial jitter
#'    applied to a sampling point in the Cartesian x- (x-axis) and y-coordinates (y-axis).
#'
#' @param boundary (Optional) An object of class SpatialPolygons (see sp::SpatialPolygons()) with
#' the outer and inner limits of the spatial sampling domain (see `candi`). These SpatialPolygons
#' help depict the spatial distribution of the (starting and current) sample configuration inside
#' the spatial sampling domain. The outer limits of `candi` serve as a rough `boundary` when
#' `plotit = TRUE`, but the SpatialPolygons are missing.
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
#' energy states return along with the results.
#'
#' @keywords spatial optimize iteration
#'
#' @concept simulated annealing
