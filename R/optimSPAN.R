# Optimization of sample configurations for variogram and spatial trend 
# identification and estimation, and for spatial interpolation
# 
# Optimize a sample configuration for variogram and spatial trend 
# identification and estimation, and for spatial interpolation. An utility
# function \emph{U} is defined so that the sample points cover, extend over,
# spread over, \bold{SPAN} the feature, variogram and geographic spaces. The
# utility function is obtained aggregating four single objective functions:
# \bold{CORR}, \bold{DIST}, \bold{PPL}, and \bold{MSSD}.
# 
# @template spJitter_doc
# @template spSANN_doc
# @template ACDC_doc
# @template MOOP_doc
# 
# @return
# \code{optimSPAN()} returns a matrix: the optimized sample configuration.
# 
# @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
# @keywords spatial optimize
# @concept simulated annealing
# @importFrom pedometrics cramer
# @importFrom pedometrics is.numint
# @importFrom pedometrics cont2cat
# @importFrom SpatialTools dist2
# @export
# @examples
# MAIN FUNCTION ################################################################
optimSPAN <-
  function (
    # OPTIM
    points, candi, x.max, x.min, y.max, y.min, iterations,
    acceptance = list(initial = 0.99, cooling = iterations / 10),
    stopping = list(max.count = iterations / 10), plotit = TRUE,
    boundary, progress = TRUE, verbose = TRUE, track = TRUE, greedy = FALSE, 
    # PPL
    lags = 7, lags.type = "exponential", lags.base = 2, cutoff = NULL, 
    criterion = "distribution", distri = NULL, pairs = FALSE, 
    # CORR and DIST
    covars, strata.type = "equal.area", use.coords = FALSE, 
    # MOOP
    weights = list(CORR = 1/4, DIST = 1/4, PPL = 1/4, MSSD = 1/4),
    nadir = list(sim = NULL, seeds = NULL, user = NULL, abs = NULL),
    utopia = list(user = NULL, abs = NULL)) {
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Check other arguments
    check <- .checkPPL(lags = lags, lags.type = lags.type, pairs = pairs,
                       lags.base = lags.base, cutoff = cutoff, 
                       criterion = criterion, distri = distri, fun = "optimPPL")
    if (!is.null(check)) stop (check, call. = FALSE)
    check <- .optimACDCcheck(candi = candi, covars = covars, 
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    check <- .checkMKV(covars = covars, eqn = eqn, vgm = vgm, 
                       krige.stat = krige.stat, candi = candi)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Prepare 'covars' and create the starting sample matrix 'sm'
    eval(.prepare_acdc_covars())
    
    # Base data
    # CORR
    pcm <- .corCORR(obj = covars, covars.type = covars.type)
    scm <- .corCORR(obj = sm, covars.type = covars.type)
    # DIST
    pop_prop <- .strataACDC(n.pts = n_pts, strata.type = strata.type, 
                            covars = covars, covars.type = covars.type)
    # PPL
    cutoff <- .cutoffPPL(cutoff = cutoff, x.max = x.max, y.max = y.max)
    lags <- .lagsPPL(lags = lags, lags.type = lags.type, cutoff = cutoff, 
                     lags.base = lags.base)
    n_lags <- length(lags) - 1
    dm_ppl <- SpatialTools::dist1(conf0[, 2:3])
    ppl <- .getPPL(lags = lags, n.lags = n_lags, dist.mat = dm_ppl, 
                   pairs = pairs)
    # MSSD
    dm_mssd <- SpatialTools::dist2(candi[, 2:3], conf0[, 2:3])
    
    # Nadir and utopia points
    
    
    
    # Base data and initial energy state (energy)
    if (covars.type == "numeric") { # Numeric covariates
      pcm <- cor(covars, use = "complete.obs")
      strata <- .numStrata(n.pts = n_pts, covars = covars,
                           strata.type = strata.type)
      nadir <- .panNadir(n.pts = n_pts, n.candi = n_candi, candi = candi,
                         n.lags = n_lags, lags = lags, criterion = criterion,
                         pre.distri = pre.distri, nadir = nadir, n.cov = n_cov,
                         covars = covars, covars.type = covars.type,
                         pcm = pcm, strata = strata)
      scm <- cor(sm, use = "complete.obs")
      energy0_acdc <- .objNum(sm = sm, n.cov = n_cov, strata = strata,
                              pcm = pcm, scm = scm, nadir = nadir,
                              weights = weights, n.pts = n_pts, utopia = utopia)
      
    } else { # Factor covariates
      pcm <- pedometrics::cramer(covars)
      pop_prop <- lapply(covars, function(x) table(x) / n_candi)
      nadir <- .panNadir(n.pts = n_pts, n.candi = n_candi, candi = candi,
                         n.lags = n_lags, lags = lags, criterion = criterion,
                         pre.distri = pre.distri, nadir = nadir, n.cov = n_cov,
                         covars = covars, covars.type, pcm, pop.prop)
      scm <- pedometrics::cramer(sm)
      energy0_acdc <- .objFac(sm = sm, pop.prop = pop_prop, nadir = nadir,
                              weights = weights, pcm = pcm, scm = scm,
                              n.pts = n_pts, n.cov = n_cov, utopia = utopia)
    }
    energy0_ppl <- energy0_ppl / attr(nadir, "ppl")
    energy0_mssd <- .calcMSSDCpp(x = dm_mssd) / attr(nadir, "mssd")

    # WEIGHTED SUM METHOD
    # The sub-arguments of 'PAN$weights' must be named so that we have a
    # guarantee that the correct weights are used.
    energy0_ppl  <- energy0_ppl * PAN$weights$PPL
    energy0_acdc <- energy0_acdc * PAN$weights$ACDC
    energy0_mssd <- energy0_mssd * PAN$weights$MSSD
    energy0 <- energy0_ppl + energy0_acdc + energy0_mssd

    # OTHER SETTINGS FOR THE SIMULATED ANNEALING ALGORITHM
    MOOP <- TRUE # this is a multi-objective optimization problem
    # ppl
    old_dm_ppl   <- dm_ppl
    new_dm_ppl   <- dm_ppl
    best_dm_ppl  <- dm_ppl
    # acdc
    old_scm      <- scm
    new_scm      <- scm
    best_scm     <- scm
    old_sm       <- sm
    new_sm       <- sm
    best_sm      <- sm
    # mssd
    old_dm_mssd  <- dm_mssd
    new_dm_mssd  <- dm_mssd
    best_dm_mssd <- dm_mssd
    # other
    count <- 0
    old_energy   <- energy0
    best_energy  <- Inf
    #energies     <- accept_probs <- vector()
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()

    # BEGIN THE ITERATIONS
    for (k in 1:iterations) {

      # Plotting and jittering #################################################
      eval(.plot_and_jitter())
      ##########################################################################
      
      # UPDATE DATASETS AND COMPUTE NEW ENERGY STATE
      # ppl: distance matrix
      new_dm_ppl <- .updatePPLCpp(new_conf[, 2:3], old_dm_ppl, wp)
      ppl <- .getPointsPerLag(lags, new_dm_ppl)
      new_energy_ppl <- .objPointsPerLag(ppl, n_lags, n_pts, criterion,
                                         pre.distri)
      new_energy_ppl <- new_energy_ppl / attr(nadir, "ppl")

      # acdc: sample and correlation matrices
      if (covars.type == "numeric") { # Numeric covariates
        #new_row <- covars[new_conf[wp, 1], ]
        #new_sm[wp, ] <- new_row
        new_sm[wp, ] <- ACDC$covars[new_conf[wp, 1], ]
        new_scm <- cor(new_sm, use = "complete.obs")
        new_energy_acdc <- .objNum(sm = new_sm, n.cov = n_cov, strata = strata,
                                   pcm = pcm, scm = new_scm, nadir = nadir,
                                   weights = ACDC$weights)
      } else { # Factor covariates
        if (covars.type == "factor") {
          #new_row <- covars[new_conf[wp, 1], ]
          #new_sm[wp, ] <- new_row
          new_sm[wp, ] <- ACDC$covars[new_conf[wp, 1], ]
          new_scm <- pedometrics::cramer(new_sm)
          new_energy_acdc <- .objFac(sm = new_sm, pop.prop = pop_prop,
                                     nadir = nadir, weights = ACDC$weights,
                                     pcm = pcm, scm = new_scm, n.pts = n_pts,
                                     n.cov = n_cov)
        }
      }

      # mssd: distance matrix
      x2 <- matrix(new_conf[wp, 2:3], nrow = 1)
      new_dm_mssd <- .updateMSSDCpp(x1 = candi[, 2:3], x2 = x2,
                                    dm = old_dm_mssd, idx = wp)
      new_energy_mssd <- .calcMSSDCpp(new_dm_mssd) / attr(nadir, "mssd")

      # WEIGHTED SUM METHOD
      # The sub-arguments of 'PAN$weights' must be named so that we have a
      # guarantee that the correct weights are used.
      new_energy_ppl  <- new_energy_ppl * PAN$weights$PPL
      new_energy_acdc <- new_energy_acdc * PAN$weights$ACDC
      new_energy_mssd <- new_energy_mssd * PAN$weights$MSSD
      new_energy      <- new_energy_ppl + new_energy_acdc + new_energy_mssd

      # Evaluate the new system configuration
      random_prob <- runif(1)
      actual_prob <- acceptance[[1]] * exp(-k / acceptance[[2]])
      if (track) accept_probs[k] <- actual_prob
      if (new_energy <= old_energy) {
        old_conf    <- new_conf
        old_energy  <- new_energy
        count       <- 0
        # ppl
        old_dm_ppl  <- new_dm_ppl
        # acdc
        old_sm      <- new_sm
        old_scm     <- new_scm
        # mssd
        old_dm_mssd <- new_dm_mssd
      } else {
        if (new_energy > old_energy & random_prob <= actual_prob) {
          old_conf    <- new_conf
          old_energy  <- new_energy
          count       <- count + 1
          # ppl
          old_dm_ppl  <- new_dm_ppl
          # acdc
          old_sm      <- new_sm
          old_scm     <- new_scm
          # mssd
          old_dm_mssd <- new_dm_mssd
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... p = ",
                random_prob, "\n")
          }
        } else {
          new_energy  <- old_energy
          new_conf    <- old_conf
          count       <- count + 1
          # ppl
          new_dm_ppl  <- old_dm_ppl
          # acdc
          new_sm      <- old_sm
          new_scm     <- old_scm
          # mssd
          new_dm_mssd <- old_dm_mssd
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... stops at",
                stopping[[1]], "\n")
          }
        }
      }

      # Best energy state
      if (track) energies[k] <- new_energy
      if (new_energy < best_energy / 1.0000001) {
        best_k           <- k
        best_conf        <- new_conf
        best_energy      <- new_energy
        best_old_energy  <- old_energy
        old_conf         <- old_conf
        # ppl
        best_dm_ppl      <- new_dm_ppl
        best_old_dm_ppl  <- old_dm_ppl
        # acdc
        best_sm          <- new_sm
        best_old_sm      <- old_sm
        best_scm         <- new_scm
        best_old_scm     <- old_scm
        # mssd
        best_dm_mssd     <- new_dm_mssd
        best_old_dm_mssd <- old_dm_mssd
      }

      # Freezing parameters
      if (count == stopping[[1]]) {
        if (new_energy > best_energy * 1.000001) {
          old_conf    <- old_conf
          new_conf    <- best_conf
          old_energy  <- best_old_energy
          new_energy  <- best_energy
          count       <- 0
          # ppl
          new_dm_ppl  <- best_dm_ppl
          old_dm_ppl  <- best_old_dm_ppl
          # acdc
          new_sm      <- best_sm
          new_scm     <- best_scm
          old_sm      <- best_old_sm
          old_scm     <- best_old_scm
          # mssd
          new_dm_mssd <- best_dm_mssd
          old_dm_mssd <- best_old_dm_mssd
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
# INTERNAL FUNCTION - COMPUTE THE NADIR VALUE ##################################
.nadirACDC <-
  function (n.pts, n.cov, n.candi, nadir, candi, covars, pcm, pop.prop, 
            covars.type, lags, n.lags, pairs, distri, criterion) {
    
    # Simulate the nadir point
    if (!is.null(nadir$sim) && !is.null(nadir$seeds)) { 
      m <- paste("simulating ", nadir$sim, " nadir values...", sep = "")
      message(m)
      
      # Set variables
      nadirDIST <- vector()
      nadirCORR <- vector()
      nadirPPL <- vector()
      nadirMSSD <- vector()
      
      # Begin the simulation
      for (i in 1:nadir$sim) {
        set.seed(nadir$seeds[i])
        pts <- sample(1:n.candi, n.pts)
        sm <- covars[pts, ]
        # CORR
        scm <- .corCORR(obj = sm, covars.type = covars.type)
        nadirCORR[i] <- .objCORR(scm = scm, pcm = pcm)
        # DIST
        nadirDIST[i] <- .objDIST(sm = sm, n.pts = n.pts, n.cov = n.cov, 
                                 pop.prop = pop.prop, covars.type = covars.type)
        # PPL
        dm <- SpatialTools::dist1(pts[, 2:3])
        ppl <- .getPPL(lags = lags, n.lags = n.lags, dist.mat = dm, 
                       pairs = pairs)
        distri <- .distriPPL(n.lags = n.lags, n.pts = n.pts, distri = distri, 
                             criterion = criterion, pairs = pairs)
        nadirPPL[i] <- .objPPL(ppl = ppl, n.lags = n.lags, n.pts = n.pts,
                               criterion = criterion, distri = distri, 
                               pairs = pairs)
        # MSSD
        dm <- SpatialTools::dist2(candi[, 2:3], pts[, 2:3])
        nadirMSSD[i] <- .objMSSD(x = dm)
      }
      
      # Prepare output
      res <- list(DIST = mean(nadirDIST), CORR = mean(nadirCORR), 
                  PPL = mean(nadirPPL), MSSD = mean(nadirMSSD))
      
    } else {
      
      # User-defined nadir values
      if (!is.null(nadir$user)) { 
        res <- list(DIST = nadir$user$DIST, CORR = nadir$user$CORR,
                    PPL = nadir$user$PPL, MSSD = nadir$user$MSSD)
        
      } else {
        if (!is.null(nadir$abs)) { 
          message("sorry but the nadir point cannot be calculated")
        }
      }
    }
    return (res)
  }
