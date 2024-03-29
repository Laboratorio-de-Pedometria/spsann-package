#' __spsann__ annealing schedule
#'
#' Set the control parameters of the annealing schedule.
#'
#' @inheritParams spJitter
#' @template spJitter_doc
#' @template schedule_doc
#'
#' @param x.max,x.min,y.max,y.min Numeric value defining the minimum and maximum quantity of random
#' noise to be added to the Cartesian x- and y-coordinates. The units are the same as of the
#' Cartesian x- and y-coordinates. If missing, they are estimated from `candi` -- or `eval.grid`,
#' in case a coarser evaluation grid is used --, `x.min` and `y.min` being set to zero, and `x.max`
#' and `y.max` being set to half the maximum distance in the Cartesian x- and y-coordinates,
#' respectively.
#'
#' @param initial.acceptance (Optional) Vector with two positive numeric values between zero (0) and
#' one (1) defining the minimum and maximum initial acceptance probabilities. The initial acceptance
#' probability is the proportion of candidate spatial sample configurations that should be accepted
#' in the first Markov chain. The optimization is stopped and a warning is issued if the value is
#' not within the predefined range. Defaults to `initial.acceptance = c(0.95, 0.99)`. (Advanced
#' users only!)
#'
#' @param initial.temperature Numeric value larger than zero (0) defining the initial temperature of
#' the system. A low `initial.temperature`, combined with a low `initial.acceptance` results in the
#' optimization algorithm to behave as a greedy optimization algorithm, i.e. only better system
#' configurations are accepted. Defaults to `initial.temperature = 0.001`.
#'
#' @param temperature.decrease Numeric value between zero (0) and one (1) used as a multiplying
#' factor to decrease the temperature at the end of each Markov chain. Defaults to
#' `temperature.decrease = 0.95`. (Advanced users only!)
#'
#' @param chains Integer value defining the maximum number of Markov chain, i.e. the number of
#' cycles of jitters at which the temperature and the size of the neighbourhood should be kept
#' constant. Defaults to `chains = 500`.
#'
#' @param chain.length Integer value larger than zero (0) defining the length of each Markov chain
#' relative to the number of samples. Defaults to `chain.length = 1`, i.e. one time the number of
#' samples.
#'
#' @param stopping Integer value defining the maximum allowable number of Markov chains without
#' improvement of the objective function value. Defaults to `stopping = ceiling(chains * 0.5)`, i.e.
#' fifty percent the maximum number of chains.
#'
#' @return
#' A list with a set of control parameters of the annealing schedule.
#'
#' @references
#' Aarts, E. H. L.; Korst, J. H. M. Boltzmann machines for travelling salesman problems. _European
#' Journal of Operational Research_S, v. 39, p. 79-95, 1989.
#'
#' Černý, V. Thermodynamical approach to the travelling salesman problem: an efficient simulation
#' algorithm. _Journal of Optimization Theory and Applications_, v. 45, p. 41-51, 1985.
#'
#' Brus, D. J.; Heuvelink, G. B. M. Optimization of sample patterns for universal kriging of
#' environmental variables. _Geoderma_, v. 138, p. 86-95, 2007.
#'
#' Kirkpatrick, S.; Gelatt, C. D.; Vecchi, M. P. Optimization by simulated annealing. _Science_, v.
#' 220, p. 671-680, 1983.
#'
#' Metropolis, N.; Rosenbluth, A. W.; Rosenbluth, M. N.; Teller, A. H.; Teller, E. Equation of state
#' calculations by fast computing machines. _The Journal of Chemical Physics_, v. 21, p. 1087-1092,
#' 1953.
#'
#' van Groenigen, J.-W.; Stein, A. Constrained optimization of spatial sampling using continuous
#' simulated annealing. _Journal of Environmental Quality_. v. 27, p. 1078-1086, 1998.
#'
#' Webster, R.; Lark, R. M. _Field sampling for environmental science and management_. London:
#' Routledge, p. 200, 2013.
#'
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[spsann]{optimACDC}}, \code{\link[spsann]{optimCORR}},
#' \code{\link[spsann]{optimDIST}}, \code{\link[spsann]{optimMKV}},
#' \code{\link[spsann]{optimMSSD}}, \code{\link[spsann]{optimPPL}},
#' \code{\link[spsann]{optimSPAN}}, \code{\link[spsann]{optimUSER}}.
#' @export
#' @examples
#' #####################################################################
#' # NOTE: The settings below are unlikely to meet your needs.         #
#' #####################################################################
#' schedule <- scheduleSPSANN(initial.temperature = 100, cellsize = 30)
# FUNCTION #########################################################################################
scheduleSPSANN <-
  function(initial.temperature = 0.001, chains = 500,
           x.max, x.min = 0, y.max, y.min = 0, cellsize,
           stopping = ceiling(chains * 0.5),
           chain.length = 1,
           temperature.decrease = 0.95,
           initial.acceptance = c(0.95, 0.99)) {

    # Initial acceptance probability
    if (any(initial.acceptance > 1) || any(initial.acceptance < 0)) {
      stop("'initial.acceptance' values must be between 0 and 1")
    }
    if (initial.acceptance[1] >= initial.acceptance[2]) {
      stop("'initial.acceptance' must be a range (min, max) of values")
    }
    if (length(initial.acceptance) != 2) {
      stop("'initial.acceptance' must be a vector of length 2")
    }

    # Initial temperature
    if (initial.temperature <= 0) {
      stop("'initial.temperature' must be larger than 0")
    }
    if (temperature.decrease >= 1 || temperature.decrease <= 0)
      stop("'temperature.decrease' must be between 0 and 1")
    if (chains < 1 || !pedometrics::isNumint(chains))
      stop("'chains' must be an integer larger than 0")
    if (chain.length < 1 || !pedometrics::isNumint(chain.length))
      stop("'chain.length' must be an integer larger than 0")
    if (stopping < 1 || !pedometrics::isNumint(stopping))
      stop("'stopping' must be an integer larger than 0")
    if (missing(x.max)) x.max <- NULL
    # if (missing(x.min)) x.min <- NULL
    if (missing(y.max)) y.max <- NULL
    # if (missing(y.min)) y.min <- NULL
    if (missing(cellsize)) {
      stop("'cellsize' must be a numeric value greater than or equal to zero")
      # cellsize <- NULL
    } else {
     if (length(cellsize) == 1) {
       cellsize <- rep(cellsize, 2)
     }
    }
    # Output
    res <- list(
      initial.acceptance = initial.acceptance, initial.temperature = initial.temperature,
      temperature.decrease = temperature.decrease, chains = chains, chain.length = chain.length,
      stopping = stopping, x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min,
      cellsize = cellsize)
    return(res)
  }
