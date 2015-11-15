#' Pareto maximum and minimum values
#' 
#' Compute the Pareto maximum and minimum values of the objective functions 
#' that compose a multi-objective optimization problem (MOOP). The Pareto 
#' maximum and minimum values are then used to transform the objective 
#' functions using the \emph{upper-lower-bound approach} so that they can be 
#' aggregated into a single \emph{utility function} using the \emph{weighted 
#' sum method}.
#' 
#' @inheritParams optimACDC
#' @inheritParams spJitter
#' 
#' @param osc A list with the optimized sample configurations (OSC). Each OSC
#' must be named after the objective function with which it has been optimized.
#' For example, \code{osc = list(CORR = osc_corr, DIST = osc_dist)}.
#' 
#' @param ... Other arguments required to compute the objective function value.
#' See the help pages of the respective objective functions to see which 
#' arguments are needed.
#' 
#' @return 
#' A data frame with the Pareto maximum and minimum values.
#' 
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[spsann]{optimACDC}}, \code{\link[spsann]{SPAN}}
#' @export
#' @examples 
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' covars <- meuse.grid[, c(1, 2, 5)]
#' 
#' # CORR
#' schedule <- scheduleSPSANN(initial.temperature = 1, chains = 1)
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
#' schedule <- scheduleSPSANN(chains = 1, initial.temperature = 30)
#' set.seed(2001)
#' osc_ppl <- optimPPL(points = 100, candi = candi, schedule = schedule)
#' 
#' # MSSD
#' schedule <- scheduleSPSANN(chains = 1, initial.temperature = 5000)
#' set.seed(2001)
#' osc_mssd <- optimMSSD(points = 100, candi = candi, schedule = schedule)
#' 
#' # Pareto
#' pareto <- minmaxPareto(osc = list(DIST = osc_dist, CORR = osc_corr,
#'                                   PPL = osc_ppl, MSSD = osc_mssd),
#'                        candi = candi, covars = covars)
#' pareto
# FUNCTION - MAIN ##############################################################
minmaxPareto <-
  function (osc, candi, covars) {
    
    if (!all(names(osc) %in% c("CORR", "DIST", "MSSD", "PPL") == TRUE)) {
      idx <- which(names(osc) %in% c("CORR", "DIST", "MSSD", "PPL") == FALSE)
      message(paste("'", names(osc)[idx], "'", 
                    " not recognized as a valid name\n", sep = ""))
    }
    
    obj <- names(osc)
    # Convert numeric covariates into factor covariates
    if (pedometrics::anyFactor(covars) && !pedometrics::allFactor(covars)) {
      id <- which(!sapply(covars, is.factor))
      message(paste("converting ", length(id), 
                    " numeric covariates into factor covariates...", 
                    sep = ""))
      covars[, id] <- 
        pedometrics::stratify(x = covars[, id], n = nrow(osc[[i]]))
    }
    
    # Compute objective function values
    obj_dist <- sapply(1:length(osc), function (i) 
      objDIST(osc[[i]]@points, covars = covars, candi = candi, 
              strata.type = osc[[i]]@objective$strata.type,
              use.coords = osc[[i]]@objective$use.coords))
    
    obj_corr <- sapply(1:length(osc), function (i) 
      objCORR(osc[[i]]@points, covars = covars, candi = candi, 
              strata.type = osc[[i]]@objective$strata.type,
              use.coords = osc[[i]]@objective$use.coords))
    
    if (any(c("PPL", "MSSD") %in% names(osc))) {
      obj_ppl <- sapply(1:length(osc), function (i) 
        objPPL(osc[[i]]@points, candi = candi, 
               lags = osc[[i]]@objective$lags,
               criterion = osc[[i]]@objective$criterion, 
               pairs = osc[[i]]@objective$pairs))
      
      obj_mssd <- sapply(1:length(osc), function (i) 
        objMSSD(osc[[i]]@points, candi = candi))
      
      # Prepare output
      res <- data.frame(
        DIST = obj_dist, CORR = obj_corr, PPL = obj_ppl, MSSD = obj_mssd, 
        row.names = c("DIST", "CORR", "PPL", "MSSD"))
    } else {
      res <- data.frame(CORR = obj_corr, DIST = obj_dist,
                        row.names = c("DIST", "CORR"))
    }
    return (res)
  }
