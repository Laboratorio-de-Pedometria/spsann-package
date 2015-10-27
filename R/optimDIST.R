#' Optimization of sample configurations for spatial trend identification
#' and estimation (II)
#'
#' Optimize a sample configuration for spatial trend identification and 
#' estimation. A criterion is defined so that the sample reproduces the 
#' marginal distribution of the covariates (\bold{DIST}).
#'
#' @template spJitter_doc
#' @template spSANN_doc
#' @template ACDC_doc
#' @template MOOP_doc
#' @template DIST_doc
#' 
#' @return
#' \code{optimDIST} returns a matrix: the optimized sample configuration.
#'  
#' \code{objDIST} returns a numeric value: the energy state of the sample
#' configuration - the objective function value.
#'
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[clhs]{clhs}}
#' @aliases optimDIST objDIST
#' @import Rcpp
#' @export
#' @examples
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' covars <- meuse.grid[, 5]
#' set.seed(2001)
#' \dontrun{
#' # This example takes more than 5 seconds to run!
#' res <- optimDIST(points = 100, candi = candi, covars = covars, 
#'                  use.coords = TRUE)
#' objSPSANN(res) # 1.6505
#' objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)
#' }
#' # Random sample
#' pts <- sample(1:nrow(candi), 5)
#' pts <- cbind(pts, candi[pts, ])
#' objDIST(points = pts, candi = candi, covars = covars, use.coords = TRUE)
# MAIN FUNCTION ################################################################
optimDIST <-
  function (points, candi, iterations = 100, 
    # DIST
    covars, strata.type = "area", use.coords = FALSE,
    # SPSANN
    x.max, x.min, y.max, y.min,
    acceptance = list(initial = 0.90, cooling = iterations / 10, 
                      by = "iterations", temperature = 5, calibrate = TRUE,
                      temperature.decrease = 0.95),
    stopping = list(max.count = iterations / 10), plotit = FALSE, track = FALSE,
    boundary, progress = TRUE, verbose = FALSE, greedy = FALSE,
    # MOOP
    weights = NULL, nadir = NULL, utopia = NULL) {
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Check other arguments
    check <- .optimACDCcheck(candi = candi, covars = covars, 
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Prepare 'covars' and create the starting sample matrix 'sm'
    eval(.prepare_acdc_covars())
    
    # Base data and initial energy state (energy)
    pop_prop <- .strataACDC(n.pts = n_pts, strata.type = strata.type,
                            covars.type = covars.type, covars = covars)
    energy0 <- .objDIST(sm = sm, pop.prop = pop_prop, n.pts = n_pts,
                        n.cov = n_cov, covars.type = covars.type)
    
    # Other settings for the simulated annealing algorithm
    old_sm <- sm
    new_sm <- sm
    best_sm <- sm
    count <- 0
    old_energy <- energy0
    best_energy <- Inf
    if (progress) {
      pb <- utils::txtProgressBar(min = 1, max = iterations, style = 3) 
    }
    time0 <- proc.time()

    if (acceptance$by == "iterations") {
      actual_prob <- acceptance$initial
    } else {
      actual_temp <- acceptance$temperature
    }
    
    # Begin the iterations
    k <- 0
    repeat {
      k <- k + 1
      
      # Plotting and jittering
      eval(.plot_and_jitter())
      
      # Update sample and correlation matrices, and energy state
      new_sm[wp, ] <- covars[new_conf[wp, 1], ]
      new_energy <- .objDIST(sm = new_sm, pop.prop = pop_prop, n.pts = n_pts, 
                             n.cov = n_cov, covars.type = covars.type)
      
      # Evaluate the new system configuration
      if (greedy) {
        random_prob <- 1
      } else {
        random_prob <- stats::runif(1)
      }
      actual_prob <- acceptance[[1]] * exp(-k / acceptance[[2]])
      if (track) accept_probs[k] <- actual_prob
      if (new_energy <= old_energy) {
        old_conf <- new_conf
        old_energy <- new_energy
        count <- 0
        old_sm <- new_sm
      } else {
        if (new_energy > old_energy & random_prob <= actual_prob) {
          old_conf <- new_conf
          old_energy <- new_energy
          count <- count + 1
          old_sm <- new_sm
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... p = ",
                random_prob, "\n")
          }
        } else {
          new_energy <- old_energy
          new_conf <- old_conf
          count <- count + 1
          new_sm <- old_sm
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... stops at",
                stopping[[1]], "\n")
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
        best_sm <- new_sm
        best_old_sm <- old_sm
      }
      
      # Freezing parameters
      if (count == stopping[[1]]) {
        if (new_energy > best_energy * 1.000001) {
          old_conf <- old_conf
          new_conf <- best_conf
          old_energy <- best_old_energy
          new_energy <- best_energy
          count <- 0
          new_sm <- best_sm
          old_sm <- best_old_sm
          cat("\n", "reached maximum count with suboptimal configuration\n")
          cat("\n", "restarting with previously best configuration\n")
          cat("\n", count, "iteration(s) with no improvement... stops at",
              stopping[[1]], "\n")
        } else {
          break
        }
      }
      if (progress) utils::setTxtProgressBar(pb, k)
      
      if (acceptance$by == "iterations") {
        actual_prob <- acceptance$initial * exp(-k / acceptance$cooling)
      } else {
        actual_temp <- actual_temp * acceptance$temperature.decrease
      }
      
      if (k == iterations) { break }
    }
    
    # Prepare output
    eval(.prepare_output())
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE ############################
# Arguments:
# sm: sample matrix
# n.pts: number of points
# n.cov: number of covariates
# pop.prop: the sampling strata and population proportion, for numeric 
#           covariates, and the population proportion, for factor covariates
# covars.type: the type of covariate (numeric or factor)
.objDIST <-
  function (sm, n.pts, n.cov, pop.prop, covars.type) {
    
    if (covars.type == "numeric") {
      
      # Count the number of points per marginal sampling strata
      count <- lapply(1:n.cov, function (i) {
        graphics::hist(sm[, i], pop.prop[[1]][[i]], plot = FALSE)$counts
      })
      
      # Compute the sample proportions
      samp.prop <- lapply(1:n.cov, function(i) count[[i]] / n.pts)
      
      # Compare the sample and population proportions
      samp.prop <- sapply(1:n.cov, function (i) 
        sum(abs(samp.prop[[i]] - pop.prop[[2]][[i]])))
      
    } else { # Factor covariates
      
      # Compute the sample proportions
      samp.prop <- lapply(sm, function(x) table(x) / n.pts)
      
      # Compare the sample and population proportions
      samp.prop <- sapply(1:n.cov, function (i)
        sum(abs(samp.prop[[i]] - pop.prop[[i]])))
    }
    
    # Compute the energy value
    energy <- sum(samp.prop)
    return (energy)
  }
# CALCULATE OBJECTIVE FUNCTION VALUE ###########################################
#' @rdname optimDIST
#' @export
objDIST <-
  function (points, candi,
    # DIST
    covars, strata.type = "area", use.coords = FALSE) {
    
    # Check other arguments
    check <- .optimACDCcheck(candi = candi, covars = covars, 
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare 'covars' and create the starting sample matrix 'sm'
    eval(.prepare_acdc_covars())
    
    # Calculate the energy state
    pop_prop <- .strataACDC(n.pts = n_pts, strata.type = strata.type,
                            covars = covars, covars.type = covars.type)
    energy <- .objDIST(sm = sm, pop.prop = pop_prop, n.pts = n_pts, 
                       n.cov = n_cov, covars.type = covars.type)
    return (energy)
  }
