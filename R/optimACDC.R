#' Optimization of sample configurations for spatial trend estimation
#'
#' Optimize a sample configuration for spatial trend estimation. A criterion is 
#' defined so that the sample reproduces the association/correlation between the
#' covariates, as well as their marginal distribution (\bold{ACDC}).
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
#' \code{optimACDC} returns a matrix: the optimized sample configuration with
#' the evolution of the energy state during the optimization as an attribute.
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
  function (points, candi, covars, strata.type = "area", iterations,
            use.coords = FALSE, x.max, x.min, y.max, y.min,
            
            weights = list(CORR = 0.5, DIST = 0.5),
            nadir = list(sim = NULL, seeds = NULL, user = NULL, abs = NULL),
            utopia = list(user = NULL, abs = NULL),
            
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            track = TRUE, boundary, progress = TRUE, verbose = TRUE, 
            greedy = FALSE) {
    
    if (!is.data.frame(covars)) covars <- as.data.frame(covars)
    
    # Check spsann arguments ###################################################
    eval(.check_spsann_arguments())
    ############################################################################
    
    check <- .optimACDCcheck(candi = candi, covars = covars, 
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Set plotting options ####################################################
    eval(.plotting_options())
    ############################################################################
    
    # Prepare points and candi #################################################
    eval(.prepare_points())
    ############################################################################
    
    # Prepare for jittering ####################################################
    eval(.prepare_jittering())
    ############################################################################
    
    # Prepare 'covars' and create the starting sample matrix 'sm' ##############
    eval(.prepare_acdc_covars())
    ############################################################################
    
    # Base data and initial energy state (energy)
    pcm <- .corCORR(obj = covars, covars.type = covars.type)
    scm <- .corCORR(obj = sm, covars.type = covars.type)
    
    if (covars.type == "numeric") { # Numeric covariates
      strata <- .numStrata(n.pts = n_pts, covars = covars, 
                           strata.type = strata.type)
      nadir <- .numNadir(n.pts = n_pts, n.cov = n_cov, n.candi = n_candi, 
                         pcm = pcm, nadir = nadir, candi = candi, 
                         covars = covars, strata = strata)
      utopia <- .numUtopia(utopia = utopia)
      energy0 <- .objNum(sm = sm, n.cov = n_cov, strata = strata, pcm = pcm, 
                         scm = scm, nadir = nadir, weights = weights, 
                         n.pts = n_pts, utopia = utopia)

    } else { # Factor covariates
      if (covars.type == "factor") {
        pop_prop <- lapply(covars, function(x) table(x) / nrow(covars))
        nadir <- .facNadir(nadir = nadir, candi = candi, n.candi = n_candi,
                           n.pts = n_pts, n.cov = n_cov, covars = covars, 
                           pop.prop = pop_prop, pcm = pcm)
        utopia <- .facUtopia(utopia = utopia)
        energy0 <- .objFac(sm = sm, pop.prop = pop_prop, nadir = nadir, 
                           weights = weights, pcm = pcm, scm = scm,
                           n.pts = n_pts, n.cov = n_cov, utopia = utopia)
      }
    }

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
    #energies <- data.frame(obj = NA, CORR = NA, DIST = NA)
    #accept_probs <- vector()
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()

    # Begin the main loop
    for (k in 1:iterations) {
      
      # Plotting and jittering #################################################
      eval(.plot_and_jitter())
      ##########################################################################
      # Update sample and correlation matrices, and energy state
      if (covars.type == "numeric") { # Numeric covariates
        new_row <- covars[new_conf[wp, 1], ]
        new_sm[wp, ] <- new_row
        new_scm <- cor(new_sm, use = "complete.obs")
        new_energy <- .objNum(sm = new_sm, n.cov = n_cov, strata = strata, 
                              pcm = pcm, scm = new_scm, nadir = nadir,
                              weights = weights, n.pts = n_pts, utopia = utopia)

      } else { # Factor covariates
        if (covars.type == "factor") {
          new_row <- covars[new_conf[wp, 1], ]
          new_sm[wp, ] <- new_row
          new_scm <- pedometrics::cramer(new_sm)
          new_energy <- .objFac(sm = new_sm, pop.prop = pop_prop, scm = new_scm, 
                                nadir = nadir, weights = weights, pcm = pcm, 
                                n.pts = n_pts, n.cov = n_cov, utopia = utopia)
        }
      }
      
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

    # Prepare output ###########################################################
    eval(.prepare_output())
    ############################################################################
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
.optimACDCcheck <-
  function (candi, covars, use.coords, strata.type) {
    
    # covars
    if (ncol(covars) < 2 && use.coords == FALSE) {
      res <- paste("'covars' must have two or more columns")
      return (res)
    }
    if (nrow(candi) != nrow(covars)) {
      res <-
        paste("'candi' and 'covars' must have the same number of rows")
      return (res)
    }
        
    # strata.type
    aa <- match(strata.type, c("area", "range"))
    if (is.na(aa)) {
      res <- paste("'strata.type = ", strata.type, "' is not supported",
                   sep = "")
      return (res)
    }
    rm(aa)
    
  }
# INTERNAL FUNCTION - BREAKS FOR NUMERIC COVARIATES ############################
# Now we define the breaks and the distribution, and return it as a list
# Quantiles now honour the fact that the data are discontinuous
# NOTE: there is a problem when the number of unique values is small (3)
# TODO: build a function is pedometrics
.numStrata <-
  function (n.pts, covars, strata.type) {
    
    # equal area strata
    if (strata.type == "area") {
      n_cov <- ncol(covars)
      probs <- seq(0, 1, length.out = n.pts + 1)
      breaks <- lapply(covars, quantile, probs, na.rm = TRUE, type = 3)
      
      # ASR: This is an old implementation
      #count <- lapply(breaks, table)
      #count <- lapply(count, as.integer)
      #count <- lapply(count, function(x) {x[2] <- x[2] + x[1] - 1; x[-1L]})
      
      # Compute the proportion of points per sampling strata
      breaks <- lapply(breaks, unique)
      count <- lapply(1:n_cov, function (i)
        hist(covars[, i], breaks[[i]], plot = FALSE)$counts)
      count <- lapply(1:n_cov, function(i) count[[i]] / sum(count[[i]]))
      strata <- list(breaks, count)

    } else {
      # equal range strata
      if (strata.type == "range") {
        n_cov <- ncol(covars)
        breaks <- lapply(1:n_cov, function(i)
          seq(min(covars[, i]), max(covars[, i]), length.out = n.pts + 1))
        d <- lapply(1:n_cov, function(i)
          SpatialTools::dist2(matrix(breaks[[i]]), matrix(covars[, i])))
        d <- lapply(1:n_cov, function(i) apply(d[[i]], 1, which.min))
        breaks <- lapply(1:n_cov, function(i) breaks[[i]] <- covars[d[[i]], i])
        breaks <- lapply(breaks, unique)
        count <- lapply(1:n_cov, function (i)
          hist(covars[, i], breaks[[i]], plot = FALSE)$counts)
        
        # Compute the proportion of points per sampling strata
        count <- lapply(1:n_cov, function (i) {count[[i]] / sum(count[[i]])})
                
        # ASR: This is an old implementation to merge null strata
        #breaks <- lapply(breaks, unique)
        #count <- lapply(1:n_cov, function (i)
        #  hist(covars[, i], breaks[[i]], plot = FALSE)$counts)
        #zero <- sapply(1:n_cov, function (i) any(count[[i]] == 0))
        #zero <- which(zero)
        #for (i in zero) {
        #  wz <- which(count[[i]] == 0)
        #  breaks[[i]] <- breaks[[i]][-(wz + 1)]
        #}
        #for (i in 1:n_cov) {
        #  mini <- min(covars[, i])
        #  maxi <- max(covars[, i])
        #  count[[i]] <- diff(breaks[[i]]) / ((maxi - mini) / n.pts)
        #}
        
        strata <- list(breaks, count)
      }
    }
    return (strata)
  }
# INTERNAL FUNCTION - NADIR FOR NUMERIC COVARIATES #############################
.numNadir <-
  function (n.pts, n.cov, n.candi, pcm, nadir, candi, covars, strata) {

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
        scm <- cor(sm, use = "complete.obs")
        
        # Count the proportion of points per sampling strata
        counts <- lapply(1:n.cov, function (i)
          hist(sm[, i], strata[[1]][[i]], plot = FALSE)$counts)
        counts <- lapply(1:n.cov, function(i) counts[[i]] / n.pts)
        counts <- sapply(1:n.cov, function (i)
          sum(abs(counts[[i]] - strata[[2]][[i]])))
        
        # Compute the nadir point
        nadirDIST[i] <- sum(counts)
        nadirCORR[i] <- sum(abs(pcm - scm))
      }  
      res <- list(DIST = mean(nadirDIST), CORR = mean(nadirCORR))
            
    } else {
      
      # User-defined nadir values
      if (!is.null(nadir$user)) { 
        res <- list(DIST = nadir$user$DIST, CORR = nadir$user$CORR)
        
      } else {
        if (!is.null(nadir$abs)) { 
          message("sorry but the nadir point cannot be calculated")
        }
        # Maximum absolute value: upper bound is equal to 100
        # For a single covariate, the lower bound is equal to 0
        # For the correlation matrix, the lower bound is always equal to 0
        #     # Absolute nadir point
        #     if (pre.distri == 1) {
        #       strata_nadir <- (2 * (n.pts - 1)) * n.cov / 100
        #     } else {
        #       if (pre.distri == 0.5) {
        #         strata_nadir <- 2 * n.pts * n.cov / 100
        #       } else {
        #         strata_nadir <- rep(0, n.pts)
        #         strata_nadir[which.min(pre.distri)] <- n.pts
        #         strata_nadir <- sum(abs(strata_nadir - pre.distri))
        #         strata_nadir <- strata_nadir * n.cov / 100
        #       }
        #     }
        #     cor_nadir <- (2 * n.cov) + sum(abs(pcm - diag(n.cov)))
        #     cor_nadir <- cor_nadir / 100
      }
    }
    return (res)
  }
# INTERNAL FUNCTION - CRITERION FOR NUMERIC COVARIATES #########################
.objNum <-
  function (sm, n.cov, strata, pcm, scm, nadir, weights, n.pts, utopia) {
    
    # Count the proportion of points per sampling strata
    counts <- lapply(1:n.cov, function (i)
      hist(sm[, i], strata[[1]][[i]], plot = FALSE)$counts)
    counts <- lapply(1:n.cov, function(i) counts[[i]] / n.pts)
    counts <- sapply(1:n.cov, function (i) 
      sum(abs(counts[[i]] - strata[[2]][[i]])))
        
    # Scale the objective function values
    obj_cont <- (sum(counts) - utopia$DIST) / (nadir$DIST - utopia$DIST)
    obj_cor <- (sum(abs(pcm - scm)) - utopia$CORR) / (nadir$CORR - utopia$CORR)
      
    # Aggregate the objective function values
    obj_cont <- obj_cont * weights$DIST
    obj_cor <- obj_cor * weights$CORR
    res <- obj_cont + obj_cor
    res <- data.frame(obj = res, CORR = obj_cor, DIST = obj_cont)
    
    return (res)
  }
# INTERNAL FUNCTION - NADIR FOR FACTOR COVARIATES ##############################
.facNadir <-
  function (nadir, candi, n.candi, n.pts, n.cov, covars, pop.prop, pcm) {
    
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
        scm <- pedometrics::cramer(sm)
        samp_prop <- lapply(sm, function(x) table(x) / n.pts)
        samp_prop <- sapply(1:n.cov, function (i)
          sum(abs(samp_prop[[i]] - pop.prop[[i]])))
        nadirDIST[i] <- sum(samp_prop)
        nadirCORR[i] <- sum(abs(pcm - scm))
      }
      
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
# INTERNAL FUNCTION - CRITERION FOR FACTOR COVARIATES ##########################
.objFac <-
  function (sm, pop.prop, nadir, weights, pcm, scm, n.pts, n.cov, utopia) {
    
    samp_prop <- lapply(sm, function(x) table(x) / n.pts)
    samp_prop <- sapply(1:n.cov, function (i)
      sum(abs(samp_prop[[i]] - pop.prop[[i]])))
    
    # Scale the objective function values
    obj_cat <- (sum(samp_prop) - utopia$DIST) / (nadir$DIST - utopia$DIST)
    obj_cor <- (sum(abs(pcm - scm)) - utopia$CORR) / (nadir$CORR - utopia$CORR)
    
    # Aggregate the objective function values
    obj_cat <- obj_cat * weights$DIST
    obj_cor <- obj_cor * weights$CORR
    res <- obj_cat + obj_cor
    res <- data.frame(obj = res, CORR = obj_cor, DIST = obj_cat)
    
    return (res)
  }
# INTERNAL FUNCTION - UTOPIA POINT FOR FACTOR COVARIATES #######################
.facUtopia <-
  function (utopia) {
    
    if (!is.null(unlist(utopia$user))) {
      utopia <- list(CORR = utopia$user$CORR, DIST = utopia$user$DIST)
      
    } else {
      message("sorry but the utopia point cannot be calculated")
    }
  }
# INTERNAL FUNCTION - UTOPIA POINT FOR NUMERIC COVARIATES ######################
.numUtopia <-
  function (utopia) {
    
    if (!is.null(unlist(utopia$user))) {
      utopia <- list(CORR = utopia$user$CORR, DIST = utopia$user$DIST)
      
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
        
    if (!is.data.frame(covars)) covars <- as.data.frame(covars)
    
    # Check arguments
    check <- .optimACDCcheck(candi = candi, covars = covars, 
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points and candi #################################################
    eval(.prepare_points())
    ############################################################################
    
    # Prepare 'covars' and create the starting sample matrix 'sm' ##############
    eval(.prepare_acdc_covars())
    ############################################################################
    
    # Base data and initial energy state
    pcm <- .corCORR(obj = covars, covars.type = covars.type)
    scm <- .corCORR(obj = sm, covars.type = covars.type)
    
    if (covars.type == "numeric") { # Numeric covariates
      strata <- .numStrata(n.pts = n_pts, covars = covars, 
                           strata.type = strata.type)
      nadir <- .numNadir(n.pts = n_pts, n.cov = n_cov, n.candi = n_candi, 
                         pcm = pcm, nadir = nadir, candi = candi, 
                         covars = covars, strata = strata)
      utopia <- .numUtopia(utopia = utopia)
      energy <- .objNum(sm = sm, n.cov = n_cov, strata = strata, pcm = pcm,
                        scm = scm, nadir = nadir, weights = weights,
                        n.pts = n_pts, utopia = utopia)
      
    } else { # Factor covariates
      if (covars.type == "factor") {
        pop_prop <- lapply(covars, function(x) table(x) / nrow(covars))
        nadir <- .facNadir(nadir = nadir, candi = candi, n.candi = n_candi,
                           n.pts = n_pts, n.cov = n_cov, covars = covars, 
                           pop.prop = pop_prop, pcm = pcm)
        utopia <- .facUtopia(utopia = utopia)
        energy <- .objFac(sm = sm, pop.prop = pop_prop, nadir = nadir, 
                           weights = weights, pcm = pcm, scm = scm,
                           n.pts = n_pts, n.cov = n_cov, utopia = utopia)
      }
    }
    return (energy)
  }
