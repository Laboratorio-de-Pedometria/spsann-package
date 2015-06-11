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
# \code{optimSPAN()} returns a matrix: the optimized sample configuration with
# the evolution of the energy state during the optimization as an attribute.
# 
# @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
# @keywords spatial optimize
# @concept simulated annealing
# @importFrom fields rdist
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
    boundary, progress = TRUE, verbose = TRUE, greedy = FALSE,
    # PPL
    lags = 7, lags.type = "exponential", lags.base = 2, cutoff = NULL, 
    criterion = "distribution", distri = NULL, pairs = FALSE, 
    # CORR and DIST
    covars, strata.type = "equal.area", use.coords = FALSE, 
    # MOOP
    weights = list(CORR = 1/4, DIST = 1/4, PPL = 1/4, MSSD = 1/4),
    nadir = list(sim = NULL, seeds = NULL, user = NULL, abs = NULL),
    utopia = list(user = NULL, abs = NULL)) {
    
    # Check arguments
    if (!is.data.frame(covars)) covars <- as.data.frame(covars) 
    
    # Check spsann arguments ###################################################
    eval(.check_spsann_arguments())
    ############################################################################
    
    check <- .optimPPLcheck(lags = lags, lags.type = lags.type, pairs = pairs,
                            lags.base = lags.base, cutoff = cutoff, 
                            criterion = criterion, distri = distri)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    check <- .optimACDCcheck(candi = candi, covars = covars,
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
#     check <- .MOOPcheck(weights = weights, nadir = nadir, utopia = utopia)
#     if (!is.null(check)) stop (check, call. = FALSE)

    # Set plotting options ####################################################
    eval(.plotting_options())
    ############################################################################

    # Prepare points and candi #################################################
    eval(.prepare_points())
    ############################################################################
    
    # Prepare for jittering ####################################################
    eval(.prepare_jittering())
    ############################################################################
    
    # BASE VARIABLES AND DATASETS, NADIR POINT AND INITIAL ENERGY STATE
    
    # PPL
    if (length(lags) >= 3) {
      n_lags <- length(lags) - 1
    } else {
      n_lags <- lags
      lags <- .getLagBreaks(lags, lags.type, cutoff, lags.base)
    }
    dm_ppl <- as.matrix(dist(conf0[, 2:3], method = "euclidean"))
    ppl <- .getPointsPerLag(lags, dm_ppl)
    energy0_ppl <- .objPointsPerLag(ppl = ppl, n.lags = n_lags, n.pts = n_pts,
                                    criterion = criterion, distri = distri)

    # MSSD
    dm_mssd <- fields::rdist(candi[, 2:3], conf0[, 2:3])

    # ACDC
    # Prepare covariates (covars) and create the starting sample matrix (sm)
    covars.type <- ifelse(pedometrics::is.any.factor(covars), "factor",
                          "numeric")
    covars <- .covarsACDC(covars = covars, covars.type = covars.type, 
                          use.coords = use.coords, candi = candi, n.pts = n_pts,
                          strata.type = strata.type)
    n_cov <- ncol(covars)
    sm <- covars[points[, 1], ]
    
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
    energies     <- accept_probs <- vector()
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()

    # BEGIN THE ITERATIONS
    for (k in 1:iterations) {

      # Jitter one of the points and update x.max and y.max
      # Which point (wp) are we jittering?
      wp <- sample(c(1:n_pts), 1)
      new_conf <- spJitterFinite(old_conf, candi, x.max, x.min, y.max,
                                 y.min, wp)
      x.max <- x_max0 - (k / iterations) * (x_max0 - x.min)
      y.max <- y_max0 - (k / iterations) * (y_max0 - y.min)

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
      accept_probs[k] <- actual_prob
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
      energies[k] <- new_energy
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

      # Plotting
      #if (plotit && any(round(seq(1, iterations, 10)) == k)) {
      if (plotit && pedometrics::is.numint(k / 10)) {
        .spSANNplot(energy0, energies, k, acceptance,
                    accept_probs, boundary, new_conf[, 2:3],
                    conf0[, 2:3], y_max0, y.max, x_max0, x.max,
                    best.energy = best_energy, best.k = best_k)
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
    if (progress) close(pb)
    res <- .spSANNout(new_conf, energy0, energies, time0)
    return (res)
  }
# INTERNAL FUNCTION - CALCULATE NADIR ##########################################
.panNadir <-
  function (n.pts, n.candi, candi, n.lags, lags, criterion, pre.distri,
            nadir, n.cov, covars, covars.type, pop.prop, pcm, strata) {

    if (!is.null(nadir[[1]]) && !is.null(nadir[[2]])) { # Simulate the nadir
      m <- paste("simulating ", nadir[[1]], " nadir values...", sep = "")
      message(m)

      # SET VARIABLES
      ppl_nadir    <- vector()
      strata_nadir <- vector()
      correl_nadir <- vector()
      mssd_nadir   <- vector()

      # BEGIN THE SIMULATION
      for (i in 1:nadir[[1]]) {
        pts <- sample(1:n.candi, n.pts)
        sm <- covars[pts, ]
        pts <- candi[pts, ]

        # ppl
        dm_ppl <- as.matrix(dist(pts[, 2:3], method = "euclidean"))
        ppl <- .getPointsPerLag(lags, dm_ppl)
        if (criterion == "distribution") {
          if (is.null(pre.distri)) {
            pre.distri <- rep(n.pts, n.lags)
          }
          ppl_nadir[i] <- sum(pre.distri - ppl)
        } else {
          if (criterion == "minimum") {
            ppl_nadir[i] <- n.pts / (min(ppl) + 1)
          }
        }

        # acdc
        if (covars.type == "numeric") { # numeric covariates
          scm <- cor(sm, use = "complete.obs",  method = "pearson")
          counts <- lapply(1:n.cov, function (i)
            hist(sm[, i], strata[[1]][[i]], plot = FALSE)$counts)
          counts <- sapply(1:n.cov, function (i)
            sum(abs(counts[[i]] - strata[[2]][[i]])))
          strata_nadir[i] <- sum(counts)
          correl_nadir[i] <- sum(abs(pcm - scm))

        } else { # factor covariates
          if (covars.type == "factor") {
            scm <- pedometrics::cramer(sm)
            samp_prop <- lapply(sm, function(x) table(x) / n.pts)
            samp_prop <- sapply(1:n.cov, function (i)
              sum(abs(samp_prop[[i]] - pop.prop[[i]])))
            strata_nadir[i] <- sum(samp_prop)
            correl_nadir[i] <- sum(abs(pcm - scm))
          }
        }

        # mssd
        dm_mssd <- fields::rdist(candi[, 2:3], pts[, 2:3])
        mssd_nadir[i] <- .calcMSSDCpp(dm_mssd)
      }

      # SHOULD WE SAVE AND RETURN THE SIMULATED VALUES?
      # ASR: We compute the mean simulated value and return it as an attribute
      #      because we want to explore the simulated values in the future.
      if (nadir[[2]]) {
        res <- list(ppl = ppl_nadir, strata = strata_nadir,
                    correl = correl_nadir, mssd_nadir = mssd_nadir)
      } else {
        res <- list(ppl = "ppl_nadir", strata = "strata_nadir",
                    correl = "correl_nadir", mssd_nadir = "mssd_nadir")
      }
      a <- attributes(res)
      a$ppl_nadir <- mean(ppl_nadir) / 100
      a$strata_nadir <- mean(strata_nadir) / 100
      a$correl_nadir <- mean(correl_nadir) / 100
      a$mssd_nadir <- mean(mssd_nadir) / 100
      attributes(res) <- a

      # ASR: Other options are not available yet.
    } else {
      if (!is.null(nadir[[3]])) {
        message("sorry but you cannot set the nadir point")
      } else {
        if (!is.null(nadir[[4]])) {
          message("sorry but the nadir point cannot be calculated")
        }
      }
    }
    return (res)
  }
