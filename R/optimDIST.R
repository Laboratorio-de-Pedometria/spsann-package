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
#' schedule <- scheduleSPSANN(initial.temperature = 0.5, chains = 1)
#' set.seed(2001)
#' \dontrun{
#' # This example takes more than 5 seconds to run!
#' res <- optimDIST(points = 100, candi = candi, covars = covars,
#'                  use.coords = TRUE, plotit = TRUE, schedule = schedule)
#' objSPSANN(res) # 2.170422
#' objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)
#' }
#' # Random sample
#' pts <- sample(1:nrow(candi), 5)
#' pts <- cbind(pts, candi[pts, ])
#' objDIST(points = pts, candi = candi, covars = covars, use.coords = TRUE)
# MAIN FUNCTION ################################################################
optimDIST <-
  function (points, candi,
            # DIST
            covars, strata.type = "area", use.coords = FALSE,
            # SPSANN
            schedule = scheduleSPSANN(), plotit = FALSE, track = FALSE,
            boundary, progress = TRUE, verbose = FALSE,
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
    old_energy <- energy0
    best_energy <- Inf
    if (progress) {
      max <- n_pts * schedule$chains * schedule$chain.length
      pb <- utils::txtProgressBar(min = 1, max = max, style = 3) 
    }
    time0 <- proc.time()
    
    actual_temp <- schedule$initial.temperature
    k <- 0 # count the number of jitters
    
    # Initiate the annealing schedule
    for (i in 1:schedule$chains) {
      n_accept <- 0
      
      for (j in 1:schedule$chain.length) { # Initiate one chain
        
        for (wp in 1:n_pts) { # Initiate loop through points
          k <- k + 1
          
          # Plotting and jittering
          eval(.plot_and_jitter())
          
          # Update sample and correlation matrices, and energy state
          new_sm[wp, ] <- covars[new_conf[wp, 1], ]
          new_energy <- 
            .objDIST(sm = new_sm, pop.prop = pop_prop, n.pts = n_pts, 
                     n.cov = n_cov, covars.type = covars.type)
          
          # Evaluate the new system configuration
          accept <- min(1, exp((old_energy - new_energy) / actual_temp))
          accept <- floor(rbinom(n = 1, size = 1, prob = accept))
          if (accept) {
            old_conf <- new_conf
            old_energy <- new_energy
            old_sm <- new_sm
            n_accept <- n_accept + 1
          } else {
            new_energy <- old_energy
            new_conf <- old_conf
            new_sm <- old_sm
          }
          if (track) energies[k] <- new_energy
          
          # Record best energy state
          if (new_energy < best_energy / 1.0000001) {
            best_k <- k
            best_conf <- new_conf
            best_energy <- new_energy
            best_old_energy <- old_energy
            old_conf <- old_conf
            best_sm <- new_sm
            best_old_sm <- old_sm
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
            new_sm <- best_sm
            old_sm <- best_old_sm
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
      x.max <- x_max0 - (i / schedule$chains) * (x_max0 - x.min)
      y.max <- y_max0 - (i / schedule$chains) * (y_max0 - y.min)
      
    } # End the annealing schedule
    
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
