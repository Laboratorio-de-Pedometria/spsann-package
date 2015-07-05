#' Optimization of sample configurations using a user-defined objective function
#'
#' Optimize a sample configuration using a user-defined objective function.
#' 
#' @template spJitter_doc
#' @template spSANN_doc
#' @template MOOP_doc
#' 
#' @param fun A function defining the objective function that should be
#' used to evaluate the energy state of the system configuration at each 
#' iteration. See \sQuote{Details} for more information.
#' 
#' @param ... Other arguments passed to the objective function. See 
#' \sQuote{Details} for more information.
#'
#' @details
#' The user-defined objective function \code{fun} must be an object of class 
#' \code{\link[base]{function}} and include the argument \code{points}. The
#' argument \code{points} is defined in \code{optimUSER} as a matrix with three 
#' columns: \code{[, 1]} the identification of each sample point given by the
#' respective row indexes of \code{candi}, \code{[, 2]} the x-coordinates, and 
#' \code{[, 3]} the y-coordinates. The identification is useful to retrieve 
#' information from any data matrix used by the objective function defined by 
#' the user.
#' 
#' @return
#' \code{optimUSER} returns a matrix: the optimized sample configuration.
#'
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @aliases optimUSER
#' @export
#' @examples
#' require(sp)
#' require(SpatialTools)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' 
#' # Define the objective function - number of points per lag distance class
#' objUSER <-
#'   function (points, lags, n_lags, n_pts) {
#'     dm <- SpatialTools::dist1(points[, 2:3])
#'     ppl <- vector()
#'     for (i in 1:n_lags) {
#'       n <- which(dm > lags[i] & dm <= lags[i + 1], arr.ind = TRUE)
#'       ppl[i] <- length(unique(c(n)))
#'     }
#'     distri <- rep(n_pts, n_lags)
#'     res <- sum(distri - ppl)
#'   }
#' lags <- seq(1, 1000, length.out = 10)
#' 
#' # Run the optimization using the user-defined objective function
#' set.seed(2001)
#' timeUSER <- Sys.time()
#' resUSER <- optimUSER(points = 100, fun = objUSER, lags = lags, n_lags = 9,
#'                      n_pts = 100, candi = candi)
#' timeUSER <- Sys.time() - timeUSER
#' 
#' # Run the optimization using the respective function implemented in spsann
#' set.seed(2001)
#' timePPL <- Sys.time()
#' resPPL <- optimPPL(points = 100, candi = candi, lags = lags)
#' timePPL <- Sys.time() - timePPL
#' 
#' # Compare results
#' timeUSER
#' timePPL
#' lapply(list(resUSER, resPPL), countPPL, lags = lags, pairs = FALSE)
#' attr(resUSER, "energy.state") # 58
#' attr(resPPL, "energy.state") # 58
# FUNCTION - MAIN ##############################################################
optimUSER <-
  function (points, candi, iterations = 100,
    # USER
    fun, ...,
    # SPSANN
    x.max, x.min, y.max, y.min,
    acceptance = list(initial = 0.99, cooling = iterations / 10),
    stopping = list(max.count = iterations / 10), plotit = FALSE, track = FALSE,
    boundary, progress = TRUE, verbose = FALSE, greedy = FALSE,
    # MOOP
    weights = NULL, nadir = NULL, utopia = NULL) {
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering 
    eval(.prepare_jittering())
    
    # Initial energy state
    energy0 <- .energyUSER(fun = fun, points = old_conf, ...)
    
    # Other settings for the simulated annealing algorithm
    count <- 0
    old_energy <- energy0
    best_energy <- Inf
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()
    
    # Begin the iterations    
    for (k in 1:iterations) {
      
      # Plotting and jittering
      eval(.plot_and_jitter())
      
      # New energy state
      new_energy <- .energyUSER(fun = fun, points = new_conf, ...)
      
      # Evaluate the new system configuration
      if (greedy) {
        random_prob <- 1
      } else {
        random_prob <- runif(1)
      }
      actual_prob <- acceptance$initial * exp(-k / acceptance$cooling)
      if (track) accept_probs[k] <- actual_prob
      if (new_energy <= old_energy) {
        old_conf <- new_conf
        old_energy <- new_energy
        count <- 0
      } else {
        if (new_energy > old_energy && random_prob <= actual_prob) {
          old_conf <- new_conf
          old_energy <- new_energy
          count <- count + 1
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... p = ",
                random_prob, "\n")
          }
        } else {
          new_energy <- old_energy
          new_conf <- old_conf
          count <- count + 1
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... stops at",
                stopping$max.count, "\n")
          }
        }
      }
      
      # Best energy state
      if (track) energies[k] <- new_energy
      if (new_energy < best_energy / 1.0000001) {
        best_k <- k
        best_conf <- new_conf
        best_energy <- new_energy
        best_old_energy <- old_energy
        old_conf <- old_conf
      }
      
      # Freezing parameters
      if (count == stopping$max.count) {
        if (new_energy > best_energy * 1.000001) {
          old_conf <- old_conf
          new_conf <- best_conf
          old_energy <- best_old_energy
          new_energy <- best_energy
          count <- 0
          cat("\n", "reached maximum count with suboptimal configuration\n")
          cat("\n", "restarting with previously best configuration\n")
          cat("\n", count, "iteration(s) with no improvement... stops at",
              stopping$max.count, "\n")
        } else {
          break
        }
      }
      if (progress) setTxtProgressBar(pb, k)
    }
    
    # Prepare output
    eval(.prepare_output())
  }
# INTERNAL FUNCTION - CALCULATE DE ENERGY STATE ################################
# fun: objective function definition
# points: system configuration
.energyUSER <- 
  function (fun, points, ...) {
    if (missing(fun) || missing(points)) {
      stop ("'fun' and 'points' are mandatory arguments")
    }
    return (do.call(fun, list(points, ...)))
  }
