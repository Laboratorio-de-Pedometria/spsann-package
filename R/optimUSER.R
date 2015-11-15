#' Optimization of sample configurations using a user-defined objective function
#'
#' Optimize a sample configuration using a user-defined objective function.
#' 
#' @inheritParams spJitter
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
#'
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
#' \dontrun{
#' # This example takes more than 5 seconds to run!
#' require(sp)
#' require(SpatialTools)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' schedule <- scheduleSPSANN(chains = 1, initial.temperature = 30)
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
#'                      n_pts = 100, candi = candi, schedule = schedule)
#' timeUSER <- Sys.time() - timeUSER
#' 
#' # Run the optimization using the respective function implemented in spsann
#' set.seed(2001)
#' timePPL <- Sys.time()
#' resPPL <- optimPPL(points = 100, candi = candi, lags = lags, 
#'                    schedule = schedule)
#' timePPL <- Sys.time() - timePPL
#' 
#' # Compare results
#' timeUSER
#' timePPL
#' lapply(list(resUSER, resPPL), countPPL, candi = candi, lags = lags)
#' objSPSANN(resUSER) # 92
#' objSPSANN(resPPL) # 92
#' }
# FUNCTION - MAIN ##############################################################
optimUSER <-
  function (points, candi,
    # USER
    fun, ...,
    # SPSANN
    schedule = scheduleSPSANN(), plotit = FALSE, track = FALSE,
    boundary, progress = TRUE, verbose = FALSE,
    # MOOP
    weights = NULL, nadir = NULL, utopia = NULL) {
    
    # Objective function name
    objective <- "USER"
    
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
    old_energy <- energy0
    best_energy <- Inf
    actual_temp <- schedule$initial.temperature
    k <- 0 # count the number of jitters
    if (progress) {
      max <- n_pts * schedule$chains * schedule$chain.length
      pb <- utils::txtProgressBar(min = 1, max = max, style = 3)
    }
    time0 <- proc.time()
    
    # Initiate the annealing schedule
    for (i in 1:schedule$chains) {
      n_accept <- 0
      
      for (j in 1:schedule$chain.length) { # Initiate one chain
        
        for (wp in 1:n_pts) { # Initiate loop through points
          k <- k + 1
      
      # Plotting and jittering
      eval(.plot_and_jitter())
      
      # New energy state
      new_energy <- .energyUSER(fun = fun, points = new_conf, ...)
      
      # Evaluate the new system configuration
      accept <- min(1, exp((old_energy - new_energy) / actual_temp))
      accept <- floor(rbinom(n = 1, size = 1, prob = accept))
      if (accept) {
        old_conf <- new_conf
        old_energy <- new_energy
        n_accept <- n_accept + 1
      } else {
        new_energy <- old_energy
        new_conf <- old_conf
      }
      if (track) energies[k] <- new_energy
      
      # Record best energy state
      if (new_energy < best_energy / 1.0000001) {
        best_k <- k
        best_conf <- new_conf
        best_energy <- new_energy
        best_old_energy <- old_energy
        old_conf <- old_conf
      }
      
      
      if (progress) utils::setTxtProgressBar(pb, k)
        } # End loop through points
        
      } # End the chain
      
      # Check the proportion of accepted swaps in the first chain
      if (i == 1) {
        x <- round(n_accept / c(n_pts * schedule$chain.length), 2)
        if (x < schedule$initial.acceptance) {
          cat("\nlow temperature: only ", x," of acceptance in the 1st chain\n", 
              sep = "")
          break
        }
      }
      
      # Count the number of chains without any change in the objective function.
      # Restart with the previously best configuration if it exists.
      if (n_accept == 0) {
        no_change <- no_change + 1
        if (no_change > schedule$stopping) {
          if (new_energy > best_energy * 1.000001) {
            old_conf <- old_conf
            new_conf <- best_conf
            old_energy <- best_old_energy
            new_energy <- best_energy
            no_change <- 0
            cat("\nrestarting with previously best configuration\n")
          } else { 
            break 
          }
        }
        if (verbose) {
          cat("\n", no_change, "chain(s) with no improvement... stops at",
              schedule$stopping, "\n")
        }
      } else {
        no_change <-  0
      }
      
      # Update control parameters
      actual_temp <- actual_temp * schedule$temperature.decrease
      x.max <- x_max0 - (i / schedule$chains) * (x_max0 - x.min) + cellsize[1]
      y.max <- y_max0 - (i / schedule$chains) * (y_max0 - y.min) + cellsize[2]
      
    } # End the annealing schedule
    
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
