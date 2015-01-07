#' Optimization of spatial samples for trend estimation
#' 
#' Optimize spatial samples for trend estimaton using spatial simulated
#' annealing.
#' 
#' @template spJitter_doc
#' @template spSANN_doc
#' @param covars Data frame or matrix with the covariates in the columns.
#' @param continuous Logical informing if the covariates are of type 
#' \sQuote{continuous} or \sQuote{categorical}. Defaults to 
#' \code{continuous = TRUE}.
#' @param weights List with two components setting the weights assigned to the
#' sampling strata/classes and the correlation/association measure. The weights
#' must sum to unity. Defaults to 
#' \code{weights = list(strata = 0.5, correl = 0.5)}.
#' @param use.coords Logical for using the geographic coordinates as covariates.
#' Defaults to \code{use.coords = FALSE}.
#' @param strata.type Character setting the type of strata to be used with 
#' continuous covariates. Available options are \code{"equal.area"} and 
#' \code{"equal.range"}. Defaults to \code{strata.type = "equal.area"}. See
#' \sQuote{Details} for more information.
#' @param sim.nadir Number of random realizations to estimate the nadir point.
#' Defaults to \code{sim.nadir = 1000}. \sQuote{Details} for more information.
#' 
#' @details
#' 
#' @return A matrix (the optimized sample points) with attributes (the evolution
#' of the energy state and the running time).
#' 
#' @references
#' Minasny, B.; McBratney, A. B. A conditioned Latin hypercube method for
#' sampling in the presence of ancillary information. \emph{Computers &
#' Geosciences}, v. 32, p. 1378-1388, 2006.
#' 
#' Minasny, B.; McBratney, A. B. Conditioned Latin Hypercube Sampling for
#' calibrating soil sensor data to soil properties. Chapter 9. Viscarra Rossel,
#' R. A.; McBratney, A. B.; Minasny, B. (Eds.) \emph{Proximal Soil Sensing}.
#' Amsterdam: Springer, p. 111-119, 2010.
#' 
#' Mulder, V. L.; de Bruin, S.; Schaepman, M. E. Representing major soil
#' variability at regional scale by constrained Latin hypercube sampling of
#' remote sensing data. \emph{International Journal of Applied Earth Observation
#' and Geoinformation}, v. 21, p. 301-310, 2013.
#' 
#' Roudier, P.; Beaudette, D.; Hewitt, A. A conditioned Latin hypercube sampling
#' algorithm incorporating operational constraints. 5th Global Workshop on
#' Digital Soil Mapping. Sydney: p. 227-231, 2012.
#' 
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' 
#' with contributions of Gerard Heuvelink \email{gerard.heuvelink@@wur.nl} and
#' Dick Brus \email{dick.brus@@wur.nl}
#' 
#' @note
#' Some of the solutions used here were found in the source code of other
#' R-packages. As such, the authors of those packages (\pkg{intamapInteractive}
#' - Edzer Pebesma <\email{edzer.pebesma@@uni-muenster.de}> and Jon Skoien
#' <\email{jon.skoien@@gmail.com}>; \pkg{clhs} - Pierre Roudier
#' <\email{roudierp@@landcareresearch.co.nz}>) are entitled 
#' \sQuote{contributors} to the R-package \pkg{pedometrics}.
#' 
#' @seealso
#' \code{\link[clhs]{clhs}}
#' @keywords spatial optimize
#' @concept simulated annealing
#' @examples
#' require(sp)
#' require(rgeos)
#' data(meuse.grid)
#' candidates <- meuse.grid[, 1:2]
#' coordinates(candidates) <- ~ x + y
#' gridded(candidates) <- TRUE
#' boundary <- as(candidates, "SpatialPolygons")
#' boundary <- gUnionCascaded(boundary)
#' candidates <- coordinates(candidates)
#' candidates <- matrix(cbind(c(1:dim(candidates)[1]), candidates), ncol = 3)
#' covars <- meuse.grid[, c(1, 2, 3, 4, 5)]
#' points <- 100
#' x.max <- diff(bbox(boundary)[1, ])
#' y.min <- x.min <- 40
#' y.max <- diff(bbox(boundary)[2, ])
#' res <- spsannTrend(points, candidates, covars, x.max = x.max, 
#'                    x.min = x.min, y.max = y.max, y.min = y.min, 
#'                    boundary = boundary, sim.nadir = 1000)
#' 
# MAIN FUNCTION ################################################################
spsannCLHS <-
  function (points, candidates, covars, continuous = TRUE,
            weights = list(strata = 0.5, correl = 0.5), use.coords = FALSE,
            strata.type = "equal.area", sim.nadir = 1000,
            x.max, x.min, y.max, y.min, iterations,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            boundary, progress = TRUE, verbose = TRUE) {
    
    # Initial checks
    if (ncol(covars) < 2) stop ("'covars' must have two or more columns")
    if (ncol(candidates) != 3) stop ("'candidates' must have three columns")
    if (!is.data.frame(covars)) covars <- as.data.frame(covars)
    if (nrow(candidates) != nrow(covars))
      stop ("'candidates' and 'covars' must have the same number of rows")
    if (sum(unlist(weights)) != 1) stop ("the 'weights' must sum to 1")
    if (plotit) {
      par0 <- par()
      on.exit(suppressWarnings(par(par0)))
    }
    
    # Prepare sample points
    if (length(points) == 1 && is.numint(points)) {
      n_pts <- points
      points <- sample(1:nrow(candidates), n_pts)
      points <- candidates[points, ]
    } else {
      n_pts <- nrow(points)
    }
    conf0 <- points
    old_conf <- conf0
    
    # Prepare covariates (covars) and create the starting sample matrix (sm)
    if (use.coords) {
      if (!continuous) {
        coords <- data.frame(candidates[, 2:3])
        breaks <- .contStrata(n_pts, coords, strata.type)[[1]]
        coords <- .cont2cat(coords, breaks)
        covars <- data.frame(covars, coords)
      } else {
        covars <- data.frame(covars, candidates[, 2:3])
      }
    }
    n_cov <- ncol(covars)
    sm <- covars[points[, 1], ]
    if (n_cov == 1) sm <- data.frame(sm)
    
    # Base data and initial energy state (energy)
    if (continuous) { # Continuous covariates
      # ASR: we should compute the true population correlation matrix (pcm)
      #      and the compare it with the sample correlation matrix (scm)
      pcm <- cor(covars, use = "complete.obs")
      strata <- .contStrata(n_pts, covars, strata.type)
      nadir <- .contNadir(n_pts, pcm, sim.nadir, candidates, covars, strata)
      scm <- cor(sm, use = "complete.obs")
      energy0 <- .objCont(sm, strata, pcm, scm, nadir, weights)
      
    } else { # Categorical covariates
      pcm <- cramer(covars)
      pop_prop <- lapply(covars, function(x) table(x) / nrow(covars))
      
      nadir <- .catNadir(sim.nadir, candidates, n_pts, covars, pop_prop, pcm)
      scm <- cramer(sm)
      energy0 <- .objCat(sm, pop_prop, nadir, weights, pcm, scm, n_pts)
    }
    
    # Other settings for the simulated annealing algorithm
    old_scm <- new_scm <- best_scm <- scm
    old_sm <- new_sm <- best_sm <- sm
    count <- 0
    old_energy <- energy0
    best_energy <- Inf
    energies <- vector()
    accept_probs <- vector()
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
      
      # Update sample and correlation matrices, and energy state
      if (continuous) { # Continuous covariates
        new_row <- covars[new_conf[wp, 1], ]
        new_sm[wp, ] <- new_row
        new_scm <- cor(new_sm, use = "complete.obs")
        new_energy <- .objCont(new_sm, strata, pcm, new_scm, nadir, weights)
        
      } else { # Categorical covariates
        new_row <- covars[new_conf[wp, 1], ]
        new_sm[wp, ] <- new_row
        new_scm <- cramer(new_sm)
        new_energy <- .objCat(new_sm, pop_prop, nadir, weights, pcm, new_scm,
                              n_pts)
      }
      
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
# INTERNAL FUNCTION - BREAKS FOR CONTINUOUS COVARIATES #########################
# Now we define the breaks and the distribution, and return it as a list
# Quantiles now honour the fact that the data are discontinuous
.contStrata <-
  function (n.pts, covars, strata.type) {
    if (strata.type == "equal.area") {
      probs <- seq(0, 1, length.out = n.pts + 1)
      breaks <- lapply(covars, quantile, probs, na.rm = TRUE, type = 1)
      count <- lapply(breaks, table)
      breaks <- lapply(count, function(x) as.double(names(x)))
      count <- lapply(count, as.integer)
      count <- lapply(count, function(x) x[-1L])
      strata <- list(breaks, count)
    }
    if (strata.type == "equal.range") {
      n_cov <- ncol(covars)
      breaks <- lapply(1:n_cov, function(i) 
        seq(min(covars[, i]), max(covars[, i]), length.out = n.pts + 1))
      count <- lapply(1:n_cov, function (i)
        hist(covars[, i], breaks[[i]], plot = FALSE)$counts)
      zero <- sapply(1:n_cov, function (i) any(count[[i]] == 0))
      zero <- which(zero)
      for (i in zero) {
        wz <- which(count[[i]] == 0)
        breaks[[i]] <- breaks[[i]][-(wz + 1)]
      }
      for (i in 1:n_cov) {
        mini <- min(covars[, i])
        maxi <- max(covars[, i])
        count[[i]] <- diff(breaks[[i]]) / ((maxi - mini) / n.pts)
      }
      strata <- list(breaks, count)
    }
    return (strata)
  }
# INTERNAL FUNCTION - DISTRIBUTION FOR CONTINUOUS COVARIATES ###################
# .contDistri <-
#   function (pre.distri, n.pts) {
#     if (length(pre.distri) > 1) {
#       if (length(pre.distri) || sum(pre.distri) != n.pts) {
#         stop(paste("'pre.distri' must be of length/sum ", n.pts, sep = ""))
#       }
#     }
#     # Sample points covering the extremes of the marginal distribution
#     if (pre.distri == 0.5) {
#       pre.distri <- rep(0, n.pts)
#       pre.distri[1] <- n.pts / 2
#       pre.distri[n.pts] <- n.pts / 2
#     }
#     return (pre.distri)
#   }
# INTERNAL FUNCTION - NADIR FOR CONTINUOUS COVARIATES ##########################
# We now always simulate the nadir point
.contNadir <-
  function (n.pts, pop.cor.mat, sim.nadir, candi, covars, strata) {
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
    #         strata_nadir <- sum(abs(strata_nadir - pre.distri))*n.cov/100
    #       }
    #     }
    #     cor_nadir <- (2 * n.cov) + sum(abs(pop.cor.mat - diag(n.cov)))
    #     cor_nadir <- cor_nadir / 100
    
    # Simulate the nadir point
    message(paste("simulating ", sim.nadir, " nadir values...", sep = ""))
    strata_nadir <- vector()
    correl_nadir <- vector()
    for (i in 1:sim.nadir) {
      pts <- sample(c(1:nrow(candi)), n.pts)
      sm <- covars[pts, ]
      if (ncol(covars) == 1) sm <- data.frame(sm)
      scm <- cor(sm, use = "complete.obs")
      counts <- lapply(1:ncol(covars), function (i)
        hist(sm[, i], strata[[1]][[i]], plot = FALSE)$counts)
      counts <- sapply(1:ncol(covars), function (i)
        sum(abs(counts[[i]] - strata[[2]][[i]])))
      strata_nadir[i] <- sum(counts)
      correl_nadir[i] <- sum(abs(pop.cor.mat - scm))
    }
    res <- list(strata = strata_nadir, correl = correl_nadir)
    a <- attributes(res)
    a$strata_nadir <- mean(strata_nadir) / 100
    a$correl_nadir <- mean(correl_nadir) / 100
    attributes(res) <- a
    return (res)
  }
# INTERNAL FUNCTION - CRITERION FOR CONTINUOUS COVARIATES ######################
.objCont <-
  function (sm, strata, pcm, scm, nadir, weights) {
    counts <- lapply(1:ncol(sm), function (i)
      hist(sm[, i], strata[[1]][[i]], plot = FALSE)$counts)
    counts <- sapply(1:ncol(sm), function (i)
      sum(abs(counts[[i]] - strata[[2]][[i]])))
    obj_cont <- sum(counts) / attr(nadir, "strata")
    obj_cont <- obj_cont * weights[[1]]
    obj_cor <- sum(abs(pcm - scm)) / attr(nadir, "correl")
    obj_cor <- obj_cor * weights[[2]]
    res <- obj_cont + obj_cor
    return (res)
  }
# INTERNAL FUNCTION - NADIR FOR CATEGORICAL COVARIATES #########################
.catNadir <-
  function (sim.nadir, candi, n.pts, covars, pop.prop, pcm) {
    message("simulating nadir values...")
    strata_nadir <- vector()
    correl_nadir <- vector()
    for (i in 1:sim.nadir) {
      pts <- sample(c(1:nrow(candi)), n.pts)
      sm <- covars[pts, ]
      if (ncol(covars) == 1) {
        sm <- data.frame(sm)
      }
      scm <- cramer(sm)
      samp_prop <- lapply(sm, function(x) table(x) / n.pts)
      samp_prop <- sapply(1:ncol(covars), function (i)
        sum(abs(samp_prop[[i]] - pop.prop[[i]])))
      strata_nadir[i] <- sum(samp_prop)
      correl_nadir[i] <- sum(abs(pcm - scm))
    }
    res <- list(strata = strata_nadir, correl = correl_nadir)
    a <- attributes(res)
    a$strata_nadir <- mean(strata_nadir) / 100
    a$correl_nadir <- mean(correl_nadir) / 100
    attributes(res) <- a
    return (res) 
  }
# INTERNAL FUNCTION - CRITERION FOR CATEGORICAL COVARIATES #####################
.objCat <-
  function (sm, pop.prop, nadir, weights, pcm, scm, n.pts) {
    samp_prop <- lapply(sm, function(x) table(x) / n.pts)
    samp_prop <- sapply(1:ncol(covars), function (i)
      sum(abs(samp_prop[[i]] - pop.prop[[i]])))
    obj_cat <- sum(samp_prop) / attr(nadir, "strata")
    obj_cat <- obj_cat * weights[[1]]
    obj_cor <- sum(abs(pcm - scm)) / attr(nadir, "correl")
    obj_cor <- obj_cor * weights[[2]]
    res <- obj_cat + obj_cor
    return (res)
  }
