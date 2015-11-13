#' Pareto maximum and minimum values
#' 
#' Compute the Pareto maximum and minimum values of the objective functions 
#' that compose a multi-objective optimization problem (MOOP). The Pareto 
#' maximum and minimum values are then used to transform the objective 
#' functions using the \emph{upper-lower-bound approach} so that they can be 
#' aggregated into a single \emph{utility function} using the \emph{weighted 
#' sum method}.
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
#' schedule <- scheduleSPSANN(initial.temperature = 1, chains = 1)
#' set.seed(2001)
#' osc_corr <- optimCORR(points = 10, candi = candi, covars = covars, 
#'                       schedule = schedule)
#' set.seed(2001)
#' osc_dist <- optimDIST(points = 10, candi = candi, covars = covars,
#'                       schedule = schedule)
#' pareto <- minmaxPareto(osc = list(DIST = osc_dist, CORR = osc_corr),
#'                        candi = candi, covars = covars)
# FUNCTION - MAIN ##############################################################
minmaxPareto <-
  function (osc, ...) {
    
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
        pedometrics::stratify(x = covars[, id], n = nrow(osc$CORR))
    }
    
    # Compute objective function values
    obj_dist <- sapply(osc, objDIST, candi = candi, covars = covars)
    obj_corr <- sapply(osc, objCORR, candi = candi, covars = covars)
    if (any(c("PPL", "MSSD") %in% names(osc))) {
      obj_ppl <- sapply(osc, objPPL, candi = candi, x.max = x.max, 
                        y.max = y.max)
      obj_mssd <- sapply(osc, objMSSD, candi = candi)
      
      # Prepare output
      res <- data.frame(DIST = obj_dist, CORR = obj_corr, PPL = obj_ppl,
                        MSSD = obj_mssd)
    } else {
      res <- data.frame(CORR = obj_corr, DIST = obj_dist)
    }
    return (res)
  }
