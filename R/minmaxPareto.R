#' Pareto minimum and maximum values
#'
#' @description
#' Compute the minimum and maximum attainable values of the objective functions that compose a
#' multi-objective combinatorial optimization problem.
#'
#' @inheritParams optimACDC
#' @inheritParams spJitter
#'
#' @param osc A list of objects of class `OptimizedSampleConfiguration` (OSC). Each OSC of the list
#' must be named after the objective function with which it has been optimized. For example,
#' `osc = list(CORR = osc_corr, DIST = osc_dist)`.
#'
#' @details
#' \subsection{Multi-objective combinatorial optimization problems}{
#' A method of solving a multi-objective combinatorial optimization problem (MOCOP) is to aggregate
#' the objective functions into a single _utility function_. In __spsann__, the aggregation
#' is performed using the _weighted sum method_, which incorporates in the weights the
#' preferences of the user regarding the relative importance of each objective function.
#'
#' The weighted sum method is affected by the relative magnitude of the different function values.
#' The objective functions implemented in __spsann__ have different units and orders of magnitude.
#' The consequence is that the objective function with the largest values may have a numerical
#' dominance during the optimization. In other words, the weights may not express the true
#' preferences of the user, resulting that the meaning of the utility function becomes unclear
#' because the optimization will favour the objective function which is numerically dominant.
#'
#' A reasonable solution to avoid the numerical dominance of any objective function is to scale the
#' objective functions so that they are constrained to the same approximate range of values. Several
#' function-transformation methods can be used for this end and __spsann__ has four of them
#' available.
#'
#' The _upper-lower-bound approach_ requires the user to inform the maximum (nadir point) and
#' minimum (utopia point) absolute function values. The resulting function values will always range
#' between 0 and 1.
#'
#' The _upper-bound approach_ requires the user to inform only the nadir point, while the utopia
#' point is set to zero. The upper-bound approach for transformation aims at equalizing only the
#' upper bounds of the objective functions. The resulting function values will always be smaller
#' than or equal to 1.
#' 
#' In most cases, the absolute maximum and minimum values of an objective function cannot be
#' calculated exactly. If the user is uncomfortable with guessing the nadir and utopia points, there
#' an option for using _numerical simulations_. It consists of computing the function value for many
#' random system configurations. The mean function value obtained over multiple simulations is used
#' to set the nadir point, while the the utopia point is set to zero. This approach is similar to
#' the upper-bound approach, but the function values will have the same orders of magnitude only at
#' the starting point of the optimization. Function values larger than one are likely to occur
#' during the optimization. We recommend the user to avoid this approach whenever possible because
#' the effect of the starting configuration on the optimization as a whole usually is insignificant
#' or arbitrary.
#'
#' The _upper-lower-bound approach_ with the minimum and maximum _attainable_ values of the
#' objective functions that compose the MOCOP, also known as the _Pareto minimum and maximum
#' values_, is the most elegant solution to scale the objective functions. However, it is the most
#' time consuming. It works as follows:
#'
#' 1. Optimize a sample configuration with respect to each objective function that composes the
#'    MOCOP;
#' 2. Compute the function value of every objective function that composes the MOCOP for every
#'    optimized sample configuration;
#' 3. Record the minimum and maximum absolute function values attained for each objective function
#'    that composes the MOCOP -- these are the Pareto minimum and maximum.
#'
#' For example, consider ACDC, a MOCOP composed of two objective functions: CORR and FREQ. The
#' minimum absolute attainable value of CORR is obtained when the sample configuration is optimized
#' with respect to only CORR, i.e. when the evaluator and generator objective functions are the same
#' (see the intersection between the second line and second column in the table below). This is
#' the Pareto minimum of CORR. It follows that the maximum absolute attainable value of CORR is
#' obtained when the sample configuration is optimized with regard to only FREQ, i.e. when the
#' evaluator function is difference from the generator function (see the intersection between the
#' first row and the second column in the table below). This is the Pareto maximum of CORR. The same
#' logic applies for finding the Pareto minimum and maximum of FREQ.
#'
#' \tabular{rll}{
#' \emph{Evaluator} \tab \emph{Generator} \tab      \cr
#'                  \tab FREQ             \tab CORR \cr
#' FREQ             \tab 0.5              \tab 8.6  \cr
#' CORR             \tab 6.4              \tab 0.3  \cr
#' }
#'
#' }
#'
#' @return
#' A data frame with the Pareto minimum and maximum values.
#'
#' @references
#' Arora, J. _Introduction to optimum design_. Waltham: Academic Press, p. 896, 2011.
#'
#' Marler, R. T.; Arora, J. S. Survey of multi-objective optimization methods for engineering.
#' _Structural and Multidisciplinary Optimization_, v. 26, p. 369-395, 2004.
#'
#' Marler, R. T.; Arora, J. S. Function-transformation methods for multi-objective optimization.
#' _Engineering Optimization_, v. 37, p. 551-570, 2005.
#'
#' Marler, R. T.; Arora, J. S. The weighted sum method for multi-objective optimization: new
#' insights. _Structural and Multidisciplinary Optimization_, v. 41, p. 853-862, 2009.
#'
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso [spsann::optimACDC()], [spsann::SPAN()]
#' @export
#' @examples
#' #####################################################################
#' # NOTE: The settings below are unlikely to meet your needs.         #
#' #####################################################################
#' if (interactive() & require(sp)) {
#'   # This example takes more than 5 seconds to run
#'   data(meuse.grid, package = "sp")
#'   # General (greedy) cooling schedule
#'   schedule <- scheduleSPSANN(
#'     initial.acceptance = c(0.01, 0.99), chains = 100,
#'     x.max = 1540, y.max = 2060, x.min = 0,
#'     y.min = 0, cellsize = 40)
#'   # CORR
#'   schedule$initial.temperature <- 0.1
#'   set.seed(2001)
#'   osc_corr <- optimCORR(
#'     points = 10, candi = meuse.grid[, 1:2],
#'     covars = meuse.grid[, c(1, 2)],
#'     schedule = schedule)
#'   # FREQ
#'   schedule$initial.temperature <- 0.1
#'   set.seed(2001)
#'   osc_dist <- optimDIST(
#'     points = 10, candi = meuse.grid[, 1:2],
#'     covars = meuse.grid[, c(1, 2)],
#'     schedule = schedule)
#'   # PPL
#'   schedule$initial.temperature <- 0.1
#'   set.seed(2001)
#'   osc_ppl <- optimPPL(
#'     points = 10, candi = meuse.grid[, 1:2],
#'     schedule = schedule)
#'   # MSSD
#'   schedule$initial.temperature <- 0.1
#'   set.seed(2001)
#'   osc_mssd <- optimMSSD(
#'     points = 10, candi = meuse.grid[, 1:2],
#'     schedule = schedule)
#'   # Pareto
#'   pareto <- minmaxPareto(
#'     osc = list(
#'       DIST = osc_dist,
#'       CORR = osc_corr,
#'       PPL = osc_ppl,
#'       MSSD = osc_mssd),
#'     candi = meuse.grid[, 1:2],
#'     covars = meuse.grid[, c(1, 2)])
#'   round(pareto, 4)
#' }
# FUNCTION - MAIN ##################################################################################
minmaxPareto <-
  function(osc, candi, covars) {
    obj <- c("CORR", "DIST", "PPL", "MSSD")
    # Check function arguments
    if (!all(names(osc) %in% obj == TRUE)) {
      idx <- which(names(osc) %in% obj == FALSE)
      message(paste0("'", names(osc)[idx], "'", " not recognized as a valid name\n"))
    } else {
      idx <- match(names(osc), obj)
      osc <- osc[idx]
    }
    # Convert numeric covariates into factor covariates
    if (pedometrics::anyFactor(covars) && !pedometrics::allFactor(covars)) {
      id <- which(!sapply(covars, is.factor))
      message(paste0("converting ", length(id), " numeric covariates into factor covariates..."))
      covars[, id] <- pedometrics::stratify(x = covars[, id], n = nrow(osc[[1]][["points"]]))
    }
    # Get objective function parameters
    strata_type <- osc$CORR$objective$strata.type
    use_coords <- osc$CORR$objective$use.coords
    # Compute objective function values
    obj_corr <- sapply(seq_along(osc), function(i) {
      objCORR(
        osc[[i]]$points, covars = covars, candi = candi, strata.type = strata_type,
        use.coords = use_coords)
    })
    obj_dist <- sapply(seq_along(osc), function(i) {
      objDIST(
        osc[[i]]$points, covars = covars, candi = candi, strata.type = strata_type,
        use.coords = use_coords)
    })
    if (all(c("PPL", "MSSD") %in% names(osc))) {
      # Compute objective function values
      # PPL
      lags <- osc$PPL$objective$lags
      criterion <- osc$PPL$objective$criterion
      pairs <- osc$PPL$objective$pairs
      obj_ppl <- sapply(seq_along(osc), function(i) {
        objPPL(
          points = osc[[i]]$points,
          candi = candi,
          lags = lags,
          criterion = criterion,
          pairs = pairs)
      })
      # MSSD
      obj_mssd <- sapply(seq_along(osc), function(i) {
        objMSSD(osc[[i]]$points, candi = candi)
      })
      # Prepare output
      res <- data.frame(
        CORR = obj_corr, DIST = obj_dist, PPL = obj_ppl, MSSD = obj_mssd,
        row.names = c("CORR", "DIST", "PPL", "MSSD"))
    } else {
      res <- data.frame(CORR = obj_corr, DIST = obj_dist, row.names = c("CORR", "DIST"))
    }
    return(res)
  }
