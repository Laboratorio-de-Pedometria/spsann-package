.panNadir <-
  function (n.pts, n.candi, candi, n.lags, lags, ppl.criterion, pre.distri,
            nadir, n.cov, covars, covar.type, pop.prop, pcm, strata) {
    
    if (!is.null(nadir[[1]])) { # Simulate the nadir point
      m <- paste("simulating ", nadir[[1]], " nadir values...", sep = "")
      message(m)
      
      # ppl
      ppl_nadir <- vector()
      for (i in 1:nadir[[1]]) {
        pts <- sample(1:n.candi, n.pts)
        pts <- candi[pts, ]
        dm <- as.matrix(dist(pts[, 2:3], method = "euclidean"))
        ppl <- .getPointsPerLag(lags, dm)
        if (ppl.criterion == "distribution") {
          if (is.null(pre.distri)) {
            pre.distri <- rep(n.pts, n.lags)
          }
          ppl_nadir[i] <- sum(pre.distri - ppl)
        } else {
          if (ppl.criterion == "minimum") {
            ppl_nadir[i] <- n.pts / (min(ppl) + 1)
          }
        }
      }
      
      # acdc
      if (covar.type == "numeric") { # numeric covariates
        strata_nadir <- vector()
        correl_nadir <- vector()
        for (i in 1:nadir[[1]]) {
          pts <- sample(c(1:n.candi), n.pts)
          sm <- covars[pts, ]
          scm <- cor(sm, use = "complete.obs")
          counts <- lapply(1:n.cov, function (i)
            hist(sm[, i], strata[[1]][[i]], plot = FALSE)$counts)
          counts <- sapply(1:n.cov, function (i)
            sum(abs(counts[[i]] - strata[[2]][[i]])))
          strata_nadir[i] <- sum(counts)
          correl_nadir[i] <- sum(abs(pop.cor.mat - scm))
        }
      } else { # factor covariates
        if (covar.type == "factor") {
          strata_nadir <- vector()
          correl_nadir <- vector()
          for (i in 1:nadir[[1]]) {
            pts <- sample(c(1:n.candi), n.pts)
            sm <- covars[pts, ]
            scm <- pedometrics::cramer(sm)
            samp_prop <- lapply(sm, function(x) table(x) / n.pts)
            samp_prop <- sapply(1:n.cov, function (i)
              sum(abs(samp_prop[[i]] - pop.prop[[i]])))
            strata_nadir[i] <- sum(samp_prop)
            correl_nadir[i] <- sum(abs(pcm - scm))
          }
        }
      }
      
      # mssd
      mssd_nadir <- vector()
      for (i in 1:nadir[[1]]) {
        pts <- sample(c(1:n.candi), n.pts)
        pts <- candi[pts, ]
        dm_mssd <- fields::rdist(candi[, 2:3], pts[, 2:3])
        mssd_nadir[i] <- .calcMSSDCpp(dm_mssd)
      }
      res <- list(ppl = ppl_nadir, strata = strata_nadir, 
                  correl = correl_nadir, mssd_nadir = mssd_nadir)
      a <- attributes(res)
      a$ppl_nadir <- mean(ppl_nadir) / 100
      a$strata_nadir <- mean(strata_nadir) / 100
      a$correl_nadir <- mean(correl_nadir) / 100
      a$mssd_nadir <- mean(mssd_nadir) / 100
      attributes(res) <- a
    } else {
      if (!is.null(nadir[[2]])) {
        message("sorry but you cannot set the nadir point")
      } else {
        message("sorry but the nadir point cannot be calculated")
      }
    }
    return (res)
  }
# MAIN FUNCTION ################################################################
optimPAN <-
  function (points, candidates, x.max, x.min, y.max, y.min, iterations = 10000,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10),
            # PPL arguments
            lags = 7, lags.type = "exponential", lags.base = 2, cutoff = NULL,
            ppl.criterion = "distribution", pre.distri = NULL,
            # ACDC arguments
            covars, covar.type = "numeric", strata.type = "equal.area", 
            use.coords = FALSE, weights.ACDC = list(strata = 0.5, correl = 0.5),
            # PAN arguments
            weights.PAN = list(PPL = 1/3, ACDC = 1/3, MSSD = 1/3),
            nadir = list(sim = 1000, user = NULL, abs = NULL), 
            # Other
            plotit = TRUE, boundary, progress = TRUE, verbose = TRUE) {
    
    # Check arguments
    if (!is.data.frame(covars)) covars <- as.data.frame(covars)
    check <- .spSANNcheck(points, candidates, x.max, x.min, y.max, y.min,
                          iterations, acceptance, stopping, plotit, boundary,
                          progress, verbose)
    if (!is.null(check)) stop (check, call. = FALSE)
    check <- .optimPPLcheck(lags, lags.type, lags.base, cutoff, ppl.criterion, 
                            pre.distri)
    if (!is.null(check)) stop (check, call. = FALSE)    
    check <- .optimACDCcheck(candidates, covars, covar.type, weights.ACDC, 
                             use.coords, strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    check <- .optimPANcheck(weights.PAN, nadir)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    if (plotit) {
      par0 <- par()
      on.exit(suppressWarnings(par(par0)))
    }
    
    # PREPARE SAMPLE POINTS
    n_candi <- nrow(candidates)
    if (length(points) == 1 && pedometrics::is.numint(points)) {
      n_pts <- points
      points <- sample(1:n_candi, n_pts)
      points <- candidates[points, ]
    } else {
      n_pts <- nrow(points)
    }
    conf0 <- points
    old_conf <- conf0
    
    # CALCULATE BASE DATASETS, NADIR POINT AND INITIAL ENERGY STATE
    # ppl
    if (length(lags) >= 3) {
      n_lags <- length(lags) - 1
    } else {
      n_lags <- lags
      lags <- .getLagBreaks(lags, lags.type, cutoff, lags.base) 
    }
    dm_ppl <- as.matrix(dist(conf0[, 2:3], method = "euclidean"))
    ppl <- .getPointsPerLag(lags, dm_ppl)
    # mssd
    dm_mssd <- fields::rdist(candidates[, 2:3], conf0[, 2:3])
    # acdc
    if (use.coords) {
      if (!continuous) {
        coords <- data.frame(candidates[, 2:3])
        breaks <- .contStrata(n_pts, coords, strata.type)[[1]]
        coords <- cont2cat(coords, breaks)
        covars <- data.frame(covars, coords)
      } else {
        covars <- data.frame(covars, candidates[, 2:3])
      }
    }
    n_cov <- ncol(covars)
    sm <- covars[points[, 1], ]
    if (continuous) { # Continuous covariates
      pcm <- cor(covars, use = "complete.obs")
      strata <- .contStrata(n_pts, covars, strata.type)
      
      nadir <- .panNadir(n.pts = n_pts, n.candi = n_candi, candi = candi, 
                         n.lags = n_lags, lags = lags, 
                         ppl.criterion = ppl.criterion, pre.distri = pre.distri,
                         nadir = nadir, n.cov = n_cov, covars = covars, 
                         covar.type = covar.type, pcm = pcm, strata = strata)
      
      scm <- cor(sm, use = "complete.obs")
      energy0
      
    } else { # Categorical covariates
      pcm <- pedometrics::cramer(covars)
      pop_prop <- lapply(covars, function(x) table(x) / nrow(covars))
      
      nadir <- .panNadir(n.pts = n_pts, n.candi = n_candi, candi = candi, 
                         n.lags = n_lags, lags = lags, 
                         ppl.criterion = ppl.criterion, pre.distri = pre.distri,
                         nadir = nadir, n.cov = n_cov, covars = covars, 
                         covar.type = covar.type, pcm = pcm, 
                         pop.prop = pop.prop)
      
      scm <- pedometrics::cramer(sm)
      energy0
    }
    
    # OTHER SETTINGS FOR THE SIMULATED ANNEALING ALGORITHM
    # ppl
    old_dm_ppl <- new_dm_ppl <- best_dm_ppl <- dm_ppl
    # acdc
    old_scm <- new_scm <- best_scm <- scm
    old_sm <- new_sm <- best_sm <- sm
    # mssd
    old_dm_mssd <- new_dm_mssd <- best_dm_mssd <- dm_mssd
    # other
    count <- 0
    old_energy <- energy0
    best_energy <- Inf
    energies <- accept_probs <- vector()
    x_max0 <- x.max
    y_max0 <- y.max
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()
    
    # Begin the main loop
    for (k in 1:iterations) {
      
      # Jitter one of the points and update x.max and y.max
      # Which point (wp)?
      wp <- sample(c(1:n_pts), 1)
      new_conf <- spJitterFinite(old_conf, candidates, x.max, x.min, y.max,
                                 y.min, wp)
      x.max <- x_max0 - (k / iterations) * (x_max0 - x.min)
      y.max <- y_max0 - (k / iterations) * (y_max0 - y.min)
      
      # UPDATE DATASETS
      # ppl: distance matrix
      new_dm_ppl <- .updatePPLCpp(new_conf[, 2:3], old_dm_ppl, wp)
      ppl <- .getPointsPerLag(lags, new_dm_ppl)
      new_energy_ppl <- .objPointsPerLag(ppl, n_lags, n_pts, criterion, 
                                         pre.distri)
      # acdc: sample and correlation matrices
      if (continuous) { # Continuous covariates
        new_row <- covars[new_conf[wp, 1], ]
        new_sm[wp, ] <- new_row
        new_scm <- cor(new_sm, use = "complete.obs")
        new_energy_acdc <- .objCont(new_sm, strata, pcm, new_scm, nadir, 
                                    weights)
      } else { # Categorical covariates
        new_row <- covars[new_conf[wp, 1], ]
        new_sm[wp, ] <- new_row
        new_scm <- pedometrics::cramer(new_sm)
        new_energy_acdc <- .objCat(new_sm, pop_prop, nadir, weights, pcm,
                                   new_scm, n_pts, n_cov)
      }
      # mssd: distance matrix
      x2 <- matrix(new_conf[wp, ], nrow = 1)
      new_dm_mssd <- .updateMSSDCpp(candidates, x2, old_dm_mssd, wp)
      new_energy_mssd <- .calcMSSDCpp(new_dm_mssd)
      
      # Evaluate the new system configuration
      random_prob <- runif(1)
      actual_prob <- acceptance[[1]] * exp(-k / acceptance[[2]])
      accept_probs[k] <- actual_prob
      if (new_energy <= old_energy) {
        old_conf <- new_conf
        old_energy <- new_energy
        old_sm <- new_sm
        old_scm <- new_scm
        count <- 0
      } else {
        if (new_energy > old_energy & random_prob <= actual_prob) {
          old_conf <- new_conf
          old_energy <- new_energy
          old_sm <- new_sm
          old_scm <- new_scm
          count <- count + 1
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... p = ", 
                random_prob, "\n")
          }
        } else {
          new_energy <- old_energy
          new_conf <- old_conf
          new_sm <- old_sm
          new_scm <- old_scm
          count <- count + 1
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... stops at",
                stopping[[1]], "\n")
          }
        }
      }
      # Best energy state
      energies[k] <- new_energy
      if (new_energy < best_energy / 1.0000001) {
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
      # Plotting
      if (plotit && any(round(seq(1, iterations, 10)) == k)) {
        .spSANNplot(energy0, energies, k, acceptance, 
                    accept_probs, boundary, new_conf[, 2:3],
                    conf0[, 2:3], y_max0, y.max, x_max0, x.max)
      }
      # Freezing parameters
      if (count == stopping[[1]]) {
        if (new_energy > best_energy * 1.000001) {
          old_conf <- old_conf
          new_conf <- best_conf
          new_sm <- best_sm
          new_scm <- best_scm
          old_energy <- best_old_energy
          new_energy <- best_energy
          old_sm <- best_old_sm
          old_scm <- best_old_scm
          count <- 0
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

# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
.optimPANcheck <-
  function (weights.PAN, nadir) {
    
    # weights.PAN
    if (!is.list(weights.PAN) || length(weights.PAN) != 3) {
      res <- paste("'weights.PAN' must be a list with three sub-arguments")
      return (res)
    }
    if (sum(unlist(weights.PAN)) != 1) {
      res <- paste("the 'weights.PAN' must sum to 1")
      return (res)
    }
    
    # nadir
    nadir <- list(sim = 1000, user = NULL, abs = NULL)
    if (!is.list(nadir) || missing(weights.PAN)) {
      res <- paste("'nadir' must be a list with three sub-arguments")
      return (res)
    }
    n <- which(sapply(nadir, is.null) == TRUE)
    if (length(nadir) == 2 && length(n) < 1) {
      res <- paste("'nadir' must have only one non-null sub-argument")
      return (res)
    }
    if (length(nadir) == 3 && length(n) < 2) {
      res <- paste("'nadir' must have only one non-null sub-argument")
      return (res)
    }
  }
