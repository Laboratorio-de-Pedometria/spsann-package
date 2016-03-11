#' Pareto minimum and maximum values
#' 
#' Compute the minimum and maximum attainable values of the objective functions that compose a multi-objective 
#' combinatorial optimization problem.
#' 
#' @inheritParams optimACDC
#' @inheritParams spJitter
#' 
#' @param osc A list of objects of class \code{OptimizedSampleConfiguration} (OSC). Each OSC of the list must
#' be named after the objective function with which it has been optimized. For example, 
#' \code{osc = list(CORR = osc_corr, DIST = osc_dist)}.
#' 
# @param ... Other arguments required to compute the objective function value. See the help pages of the 
# respective objective functions to see which arguments are needed.
#' 
#' @details 
#' The minimum and maximum attainable values of the objective functions that compose a multi-objective 
#' combinatorial optimization problem (MOCOP) are also known as the Pareto minimum and maximum values. 
#' The Pareto minimum and maximum values are used to scale the objective functions to the same approximate 
#' range of values using the \emph{upper-lower-bound approach} so that they can be correctly aggregated into a 
#' single \emph{utility function} using the \emph{weighted sum method}. Scalling the objective functions to 
#' the same approximate range of values is necessary to avoid the numerical dominance of those objective 
#' functions whose values have the largest magnitude among all objective functions. Eliminating the chance of
#' numerical dominance guarantees that the weights will play the wanted role during the optimization.
#' 
#' @return 
#' A data frame with the Pareto minimum and maximum values.
#' 
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[spsann]{optimACDC}}, \code{\link[spsann]{SPAN}}
#' @export
#' @examples 
#' \dontrun{
#' # This example takes more than 5 seconds
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' covars <- meuse.grid[, c(1, 2)]
#' 
#' # CORR
#' schedule <- scheduleSPSANN(initial.acceptance = 0.1, chains = 1, 
#'                            x.max = 1540, y.max = 2060, x.min = 0,
#'                            y.min = 0, cellsize = 40)
#' set.seed(2001)
#' osc_corr <- optimCORR(points = 10, candi = candi, covars = covars, 
#'                       schedule = schedule)
#' 
#' # DIST
#' set.seed(2001)
#' osc_dist <- optimDIST(points = 10, candi = candi, covars = covars,
#'                       schedule = schedule)
#' 
#' # PPL
#' set.seed(2001)
#' osc_ppl <- optimPPL(points = 10, candi = candi, schedule = schedule)
#' 
#' # MSSD
#' set.seed(2001)
#' osc_mssd <- optimMSSD(points = 10, candi = candi, schedule = schedule)
#' 
#' # Pareto
#' pareto <- minmaxPareto(osc = list(DIST = osc_dist, CORR = osc_corr,
#'                                   PPL = osc_ppl, MSSD = osc_mssd),
#'                        candi = candi, covars = covars)
#' pareto
#' }
# FUNCTION - MAIN ##############################################################
minmaxPareto <-
  function (osc, candi, covars) {
    
    obj <- c("CORR", "DIST", "PPL", "MSSD")
    if (!all(names(osc) %in% obj == TRUE)) {
      idx <- which(names(osc) %in% obj == FALSE)
      message(paste("'", names(osc)[idx], "'", " not recognized as a valid name\n", sep = ""))
    } else {
      idx <- match(names(osc), obj)
      osc <- osc[idx]
    }
    
    # Convert numeric covariates into factor covariates
    if (pedometrics::anyFactor(covars) && !pedometrics::allFactor(covars)) {
      id <- which(!sapply(covars, is.factor))
      message(paste("converting ", length(id), " numeric covariates into factor covariates...", sep = ""))
      covars[, id] <- pedometrics::stratify(x = covars[, id], n = nrow(osc[[1]][["points"]]))
    }
    
    # Get objective function parameters
    strata.type <- osc$CORR$objective$strata.type
    use.coords <- osc$CORR$objective$use.coords
    
    # Compute objective function values
    obj_corr <- sapply(1:length(osc), function (i) objCORR(
      osc[[i]]$points, covars = covars, candi = candi, strata.type = strata.type, use.coords = use.coords))
    obj_dist <- sapply(1:length(osc), function (i) objDIST(
      osc[[i]]$points, covars = covars, candi = candi, strata.type = strata.type, use.coords = use.coords))
    
    if (all(c("PPL", "MSSD") %in% names(osc))) {
      
      # Get objective function parameters
      lags <- osc$PPL$objective$lags
      criterion <- osc$PPL$objective$criterion
      pairs <- osc$PPL$objective$pairs
      # x.max <- osc$PPL@spsann$jitter$x[2]
      # x.min <- osc$PPL@spsann$jitter$x[1]
      # y.max <- osc$PPL@spsann$jitter$y[2]
      # y.min <- osc$PPL@spsann$jitter$y[1]
      
      # Compute objective function values
      obj_ppl <- sapply(1:length(osc), function (i) 
        objPPL(osc[[i]]$points, candi = candi, lags = lags, criterion = criterion, pairs = pairs))
               # ,x.max = x.max, y.max = y.max, x.min = x.min, y.min = y.min
      obj_mssd <- sapply(1:length(osc), function (i) objMSSD(osc[[i]]$points, candi = candi))
      
      # Prepare output
      res <- data.frame(
        CORR = obj_corr, DIST = obj_dist, PPL = obj_ppl, MSSD = obj_mssd, 
        row.names = c("CORR", "DIST", "PPL", "MSSD"))
    } else {
      res <- data.frame(CORR = obj_corr, DIST = obj_dist, row.names = c("CORR", "DIST"))
    }
    return (res)
  }
