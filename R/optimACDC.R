#' Optimization of sample configurations for spatial trend identification and
#' estimation
#'
#' Optimize a sample configuration for spatial trend identification and 
#' estimation. A criterion is defined so that the sample reproduces the 
#' bivariate association/correlation between the covariates, as well as their 
#' marginal distribution (\bold{ACDC}).
#'
#' @template spJitter_doc
#' @template spSANN_doc
#' @template ACDC_doc
#' @template MOOP_doc
#' 
#' @details
#' See \code{optimDIST} and \code{optimCORR}.
#' 
#' @return
#' \code{optimACDC} returns a matrix: the optimized sample configuration.
#' 
#' \code{objACDC} returns a numeric value: the energy state of the sample
#' configuration - the objective function value.
#'
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[clhs]{clhs}}, \code{\link[pedometrics]{cramer}}
#' @importFrom pedometrics cramer
#' @importFrom pedometrics is.numint
#' @importFrom pedometrics cont2cat
#' @importFrom pedometrics is.all.factor
#' @importFrom pedometrics is.any.factor
#' @importFrom SpatialTools dist2
#' @export
#' @examples
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' nadir <- list(sim = 10, seeds = 1:10)
#' utopia <- list(user = list(DIST = 0, CORR = 0))
#' covars <- meuse.grid[, 5]
#' set.seed(2001)
#' res <- optimACDC(points = 100, candi = candi, covars = covars, nadir = nadir,
#'                  use.coords = TRUE, iterations = 100, utopia = utopia, 
#'                  verbose = FALSE)
#' tail(attr(res, "energy")$obj, 1) # 0.5272031
#' objACDC(points = res, candi = candi, covars = covars, use.coords = TRUE, 
#'         nadir = nadir, utopia = utopia)
#' # MARGINAL DISTRIBUTION
#' par(mfrow = c(3, 3))
#' # Covariates
#' i <- sample(1:nrow(candi), 100)
#' hist(candi[, 1], breaks = 10)
#' hist(candi[, 2], breaks = 10)
#' hist(covars, breaks = 10)
#' # Optimized sample
#' hist(candi[res[, 1], 1], breaks = 10)
#' hist(candi[res[, 1], 2], breaks = 10)
#' hist(covars[res[, 1]], breaks = 10)
#' # Random sample
#' hist(candi[i, 1], breaks = 10)
#' hist(candi[i, 2], breaks = 10)
#' hist(covars[i], breaks = 10)
#' 
#' # LINEAR CORRELATION
#' # Covariates
#' cor(cbind(candi[, 1], candi[, 2], covars))
#' # Optimized sample
#' cor(cbind(candi[res[, 1], 1], candi[res[, 1], 2], covars[res[, 1]]))
#' # Random sample
#' cor(cbind(candi[i, 1], candi[i, 2], covars[i]))
# MAIN FUNCTION ################################################################
optimACDC <-
  function (
    # DIST and/or CORR
    covars, strata.type = "area", use.coords = FALSE, 
    # MOOP
    weights = list(CORR = 0.5, DIST = 0.5),
    nadir = list(sim = NULL, seeds = NULL, user = NULL, abs = NULL),
    utopia = list(user = NULL, abs = NULL),
    # SPSANN
    points, candi, iterations, x.max, x.min, y.max, y.min,
    acceptance = list(initial = 0.99, cooling = iterations / 10),
    stopping = list(max.count = iterations / 10), plotit = TRUE,
    track = TRUE, boundary, progress = TRUE, verbose = TRUE, greedy = FALSE) {
    
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
    
    # Base data and initial energy state
    pcm <- .corCORR(obj = covars, covars.type = covars.type)
    scm <- .corCORR(obj = sm, covars.type = covars.type)
    pop_prop <- .strataACDC(n.pts = n_pts, strata.type = strata.type, 
                            covars = covars, covars.type = covars.type)
    nadir <- .nadirACDC(n.pts = n_pts, n.cov = n_cov, n.candi = n_candi, 
                        pcm = pcm, nadir = nadir, covars.type = covars.type,
                        covars = covars, pop.prop = pop_prop, candi = candi)
    utopia <- .utopiaACDC(utopia = utopia)
    energy0 <- .objACDC(sm = sm, n.cov = n_cov, pop.prop = pop_prop, pcm = pcm, 
                        scm = scm, nadir = nadir, weights = weights, 
                        n.pts = n_pts, utopia = utopia, 
                        covars.type = covars.type)

    # Other settings for the simulated annealing algorithm
    old_scm <- scm
    new_scm <- scm
    best_scm <- scm
    old_sm <- sm
    new_sm <- sm
    best_sm <- sm
    count <- 0
    old_energy <- energy0
    best_energy <- data.frame(obj = Inf, CORR = Inf, DIST = Inf)
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()

    # Begin the main loop
    for (k in 1:iterations) {
      
      # Plotting and jittering
      eval(.plot_and_jitter())
      
      # Update sample and correlation matrices, and energy state
      new_sm[wp, ] <- covars[new_conf[wp, 1], ]
      new_scm <- .corCORR(obj = new_sm, covars.type = covars.type)
      new_energy <- .objACDC(sm = new_sm, pop.prop = pop_prop, scm = new_scm,
                             nadir = nadir, weights = weights, pcm = pcm, 
                             n.pts = n_pts, n.cov = n_cov, utopia = utopia,
                             covars.type = covars.type)
      
      # Evaluate the new system configuration
      if (greedy) {
        random_prob <- 1
      } else {
        random_prob <- runif(1)
      }
      actual_prob <- acceptance[[1]] * exp(-k / acceptance[[2]])
      if (track) accept_probs[k] <- actual_prob
      if (new_energy[1] <= old_energy[1]) {
        old_conf <- new_conf
        old_energy <- new_energy
        count <- 0
        old_sm <- new_sm
        old_scm <- new_scm
      } else {
        if (new_energy[1] > old_energy[1] & random_prob <= actual_prob) {
          old_conf <- new_conf
          old_energy <- new_energy
          count <- count + 1
          old_sm <- new_sm
          old_scm <- new_scm
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... p = ",
                random_prob, "\n")
          }
        } else {
          new_energy <- old_energy
          new_conf <- old_conf
          count <- count + 1
          new_sm <- old_sm
          new_scm <- old_scm
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... stops at",
                stopping[[1]], "\n")
          }
        }
      }
      
      # Best energy state
      if (track) energies[k, ] <- new_energy
      if (new_energy[1] < best_energy[1] / 1.0000001) {
        best_k <- k
        best_conf <- new_conf
        best_energy <- new_energy
        best_old_energy <- old_energy
        old_conf <- old_conf
        best_sm <- new_sm
        best_old_sm <- old_sm
        best_scm <- new_scm
        best_old_scm <- old_scm
      }
      
      # Freezing parameters
      if (count == stopping[[1]]) {
        if (new_energy[1] > best_energy[1] * 1.000001) {
          old_conf <- old_conf
          new_conf <- best_conf
          old_energy <- best_old_energy
          new_energy <- best_energy
          count <- 0
          new_sm <- best_sm
          new_scm <- best_scm
          old_sm <- best_old_sm
          old_scm <- best_old_scm
          cat("\n", "reached maximum count with suboptimal configuration\n")
          cat("\n", "restarting with previously best configuration\n")
          cat("\n", count, "iteration(s) with no improvement... stops at",
              stopping[[1]], "\n")
        } else {
          break
        }
      }
      if (progress) setTxtProgressBar(pb, k)
    }

    # Prepare output
    eval(.prepare_output())
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
# candi: candidate locations
# covars: covariates
# use.coords: should the coordinates be used
# strata.type: type of stratification of numeric covariates
.optimACDCcheck <-
  function (candi, covars, use.coords, strata.type) {
    
    # covars
    if (is.vector(covars)) {
      if (use.coords == FALSE) {
        res <- "'covars' must have two or more columns"
        return (res)
      }
      if (nrow(candi) != length(covars)) {
        res <- "'candi' and 'covars' must have the same number of rows"
        return (res)
      }
    } else {
      if (nrow(candi) != nrow(covars)) {
        res <- "'candi' and 'covars' must have the same number of rows"
        return (res)
      }
    }
    
    # strata.type
    aa <- match(strata.type, c("area", "range"))
    if (is.na(aa)) {
      res <- paste("'strata.type = ", strata.type, "' is not supported",
                   sep = "")
      return (res)
    }
  }
# INTERNAL FUNCTION - BREAKS FOR NUMERIC COVARIATES ############################
# Now we define the breaks and the distribution, and return it as a list.
# Quantiles now honour the fact that the data are discontinuous.
# NOTE: there might be a problem when the number of unique values is small (3)
.strataACDC <-
  function (n.pts, covars, strata.type, covars.type) {
    
    n_cov <- ncol(covars)
    
    if (covars.type == "factor") {
      
      # Compute the proportion of population points per marginal factor level
      res <- lapply(covars, function(x) table(x) / nrow(covars))
      
    } else { # Numeric covariates
      
      # equal area strata
      if (strata.type == "area") {
        
        # Compute the break points (discrete sample quantiles)
        probs <- seq(0, 1, length.out = n.pts + 1)
        breaks <- lapply(covars, quantile, probs, na.rm = TRUE, type = 3)
        
      } else { # equal range strata
        
        # Compute the break points
        breaks <- lapply(1:n_cov, function(i)
          seq(min(covars[, i]), max(covars[, i]), length.out = n.pts + 1))
        
        # Find and replace by the closest population value
        d <- lapply(1:n_cov, function(i)
          SpatialTools::dist2(matrix(breaks[[i]]), matrix(covars[, i])))
        d <- lapply(1:n_cov, function(i) apply(d[[i]], 1, which.min))
        breaks <- lapply(1:n_cov, function(i) breaks[[i]] <- covars[d[[i]], i])
      }
      
      # Keep only the unique break points
      breaks <- lapply(breaks, unique)
      
      # Compute the proportion of population points per marginal sampling strata
      count <- lapply(1:n_cov, function (i)
        hist(covars[, i], breaks[[i]], plot = FALSE)$counts)
      prop <- lapply(1:n_cov, FUN = function (i) {count[[i]] / sum(count[[i]])})
      
      # Output
      res <- list(breaks = breaks, prop = prop)  
    }
    
    return (res)
  }
# INTERNAL FUNCTION - COMPUTE THE NADIR VALUE ##################################
.nadirACDC <-
  function (n.pts, n.cov, n.candi, nadir, candi, covars, pcm, pop.prop, 
            covars.type) {
    
    # Simulate the nadir point
    if (!is.null(nadir$sim) && !is.null(nadir$seeds)) { 
      m <- paste("simulating ", nadir$sim, " nadir values...", sep = "")
      message(m)
      
      # Set variables
      nadirDIST <- vector()
      nadirCORR <- vector()
      
      # Begin the simulation
      for (i in 1:nadir$sim) {
        set.seed(nadir$seeds[i])
        pts <- sample(1:n.candi, n.pts)
        sm <- covars[pts, ]
        scm <- .corCORR(obj = sm, covars.type = covars.type)
        nadirDIST[i] <- .objDIST(sm = sm, n.pts = n.pts, n.cov = n.cov, 
                                 pop.prop = pop.prop, covars.type = covars.type)
        nadirCORR[i] <- .objCORR(scm = scm, pcm = pcm)
      }
      
      # Prepare output
      res <- list(DIST = mean(nadirDIST), CORR = mean(nadirCORR))
            
    } else {
      
      # User-defined nadir values
      if (!is.null(nadir$user)) { 
        res <- list(DIST = nadir$user$DIST, CORR = nadir$user$CORR)
        
      } else {
        if (!is.null(nadir$abs)) { 
          message("sorry but the nadir point cannot be calculated")
        }
      }
    }
    return (res)
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE ############################
# This function is used to calculate the criterion value of ACDC.
# It calculates, scales, weights, and aggregates the objective function values.
# Scaling is done using the upper-lower bound approach.
# Aggregation is done using the weighted sum method.
.objACDC <-
  function (sm, n.cov, nadir, weights, n.pts, utopia, pcm, scm, covars.type,
            pop.prop) {
    
    # DIST
    obj_dist <- .objDIST(sm = sm, n.pts = n.pts, n.cov = n.cov, 
                         pop.prop = pop.prop, covars.type = covars.type)
    obj_dist <- (obj_dist - utopia$DIST) / (nadir$DIST - utopia$DIST)
    obj_dist <- obj_dist * weights$DIST
    
    # CORR
    obj_cor <- .objCORR(scm = scm, pcm = pcm)
    obj_cor <- (obj_cor - utopia$CORR) / (nadir$CORR - utopia$CORR)
    obj_cor <- obj_cor * weights$CORR
    
    # Prepare output
    res <- data.frame(obj = obj_dist + obj_cor, CORR = obj_cor, DIST = obj_dist)
    return (res)
  }
# INTERNAL FUNCTION - PREPARE THE UTOPIA POINT #################################
.utopiaACDC <-
  function (utopia) {

    if (!is.null(unlist(utopia$user))) {
      list(CORR = utopia$user$CORR, DIST = utopia$user$DIST)
      
    } else {
      message("sorry but the utopia point cannot be calculated")
    }
  }
# CALCULATE OBJECTIVE FUNCTION VALUE ###########################################
#' @rdname optimACDC
#' @export
objACDC <-
  function (points, candi, covars, strata.type = "area", 
            weights = list(CORR = 0.5, DIST = 0.5), use.coords = FALSE, 
            utopia = list(user = NULL, abs = NULL),
            nadir = list(sim = NULL, seeds = NULL, user = NULL, abs = NULL)) {
    
    # Check arguments
    check <- .optimACDCcheck(candi = candi, covars = covars, 
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare 'covars' and create the starting sample matrix 'sm'
    eval(.prepare_acdc_covars())
    
    # Compute base data
    pcm <- .corCORR(obj = covars, covars.type = covars.type)
    scm <- .corCORR(obj = sm, covars.type = covars.type)
    pop_prop <- .strataACDC(n.pts = n_pts, strata.type = strata.type,
                            covars = covars, covars.type = covars.type)
    nadir <- .nadirACDC(n.pts = n_pts, n.cov = n_cov, n.candi = n_candi, 
                        pcm = pcm, nadir = nadir, candi = candi, 
                        covars = covars, pop.prop = pop_prop, 
                        covars.type = covars.type)
    utopia <- .utopiaACDC(utopia = utopia)
    
    # Compute the energy state
    energy <- .objACDC(sm = sm, pop.prop = pop_prop, covars.type = covars.type, 
                       weights = weights, pcm = pcm, scm = scm, n.pts = n_pts, 
                       n.cov = n_cov, utopia = utopia, nadir = nadir)
    return (energy)
  }
