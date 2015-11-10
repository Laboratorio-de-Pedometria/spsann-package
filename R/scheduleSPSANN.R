#' \pkg{spsann} cooling schedule
#' 
#' Set the control parameters for the cooling schedule of \pkg{spsann} 
#' functions.
#' 
#' @inheritParams spJitter
#' 
#' @param x.max,x.min,y.max,y.min Numeric value defining the minimum and 
#' maximum quantity of random noise to be added to the projected x- and 
#' y-coordinates. The minimum quantity should be equal to, at least, the 
#' minimum distance between two neighbouring candidate locations. The units 
#' are the same as of the projected x- and y-coordinates. If missing, they 
#' are estimated from \code{candi}.
#' 
#' @param initial.acceptance Numeric value between 0 and 1 defining the initial
#' acceptance probability. Defaults to \code{initial.acceptance = 0.95}.
#' 
#' @param initial.temperature Numeric value larger than 0 defining the initial
#' temperature of the system. Defaults to \code{initial.temperature = 0.001}.
#' 
#' @param temperature.decrease Numeric value between 0 and 1 defining the factor
#' by which the temperature is decreased at the end of each Markov chain. 
#' Defaults to \code{temperature.decrease = 0.95}.
#' 
#' @param chains Integer value defining the maximum number of Markov chains.
#' Defaults to \code{chains = 500}.
#' 
#' @param chain.length Integer value defining the length of each Markov chain
#' relative to the number of points. Defaults to \code{chain.length = 1}, i.e.
#' one time the number of points.
#' 
#' @param stopping Integer value defining the maximum allowable number of 
#' Markov chains without improvement of the objective function value. Defaults 
#' to \code{stopping = 10}.
#' 
#' @return
#' A list with a set of control parameters of the cooling schedule.
#' 
#' @references 
#' Aarts, E. H. L. & Korst, J. H. M. Boltzmann machines for travelling salesman 
#' problems. \emph{European Journal of Operational Research}. v. 39, p. 79-95,
#' 1989.
#' 
#' @author 
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[spsann]{optimACDC}}, \code{\link[spsann]{optimCORR}},
#' \code{\link[spsann]{optimDIST}}, \code{\link[spsann]{optimMKV}},
#' \code{\link[spsann]{optimMSSD}}, \code{\link[spsann]{optimPPL}},
#' \code{\link[spsann]{optimSPAN}}, \code{\link[spsann]{optimUSER}}.
#' @export
#' @examples
#' schedule <- scheduleSPSANN()
# FUNCTION #####################################################################
scheduleSPSANN <-
  function (initial.acceptance = 0.95, initial.temperature = 0.001,
            temperature.decrease = 0.95, chains = 500, chain.length = 1,
            stopping = 10, x.max, x.min, y.max, y.min, cellsize) {
    
    if (initial.acceptance > 1 || initial.acceptance < 0)
      stop ("'initial.acceptance' must be between 0 and 1")
    if (initial.temperature <= 0)
      stop ("'initial.temperature' must be larger than 0")
    if (temperature.decrease >= 1 || temperature.decrease <= 0)
      stop ("'temperature.decrease' must be between 0 and 1")
    if (chains < 1 || !pedometrics::isNumint(chains))
      stop ("'chains' must be an integer larger than 0")
    if (chain.length < 1 || !pedometrics::isNumint(chain.length))
      stop ("'chain.length' must be an integer larger than 0")
    if (stopping < 1 || !pedometrics::isNumint(stopping))
      stop ("'stopping' must be an integer larger than 0")
    if (missing(x.max)) x.max <- NULL
    if (missing(x.min)) x.min <- NULL
    if (missing(y.max)) y.max <- NULL
    if (missing(y.min)) y.min <- NULL
    if (missing(cellsize)) cellsize <- NULL
    
    # Output
    res <- list(initial.acceptance = initial.acceptance, 
                initial.temperature = initial.temperature,
                temperature.decrease = temperature.decrease, 
                chains = chains, chain.length = chain.length, 
                stopping = stopping, x.max = x.max, x.min = x.min, 
                y.max = y.max, y.min = y.min, cellsize = cellsize)
    return (res)
  }
