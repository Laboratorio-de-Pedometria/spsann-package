#' Optimization of sample configurations for spatial trend estimation
#'
#' Optimize a sample configuration for trend estimation. A criterion is defined
#' so that the sample reproduces the association/correlation between the 
#' covariates (\bold{CORR}).
#'
#' @template spJitter_doc
#' @template spSANN_doc
#' @template ACDC_doc
#' 
#' @details
#' The correlation between two numeric covariates is measured using the 
#' Pearson's r, a descriptive statistic that ranges from -1 to +1. 
#' This statistic is also known as the linear correlation coefficient.
#' 
#' When the set of covariates includes factor covariates, any numeric covariate 
#' is transformed into a factor covariate. The numeric covariates are 
#' categorized using the sampling strata defined using one of the two methods 
#' available (equal-area or equal-range strata) (see more details at
#' \code{optimDist}).
#' 
#' The association between two factor covariates is measured using the Cramér's 
#' v, a descriptive statistic that ranges from 0 to 1. The closer to 1 the 
#' Cramér's v is, the stronger the association between two factor covariates. 
#' The main weakness of using the Cramér's v is that, while the Pearson's r 
#' shows the degree and direction of the association between two covariates 
#' (negative or positive), the Cramér's v only measures the degree (weak or 
#' strong).
#'
#' @return
#' \code{optimCORR} returns a matrix: the optimized sample configuration with
#' the evolution of the energy state during the optimization as an attribute.
#' 
#' \code{objCORR} returns a numeric value: the energy state of the sample
#' configuration - the objective function value.
#'
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[clhs]{clhs}}, \code{\link[pedometrics]{cramer}}
#' @importFrom pedometrics cramer
#' @importFrom pedometrics is.numint
#' @importFrom pedometrics cont2cat
#' @importFrom SpatialTools dist2
#' @export
#' @examples
#' require(pedometrics)
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' covars <- meuse.grid[, 5]
#' set.seed(2001)
#' res <- optimCORR(points = 100, candi = candi, covars = covars, 
#'                  use.coords = TRUE, iterations = 100)
#' tail(attr(res, "energy"), 1) # 0.06386069
#' objCORR(points = res, candi = candi, covars = covars, use.coords = TRUE)
# MAIN FUNCTION ################################################################
optimCORR <-
  function (points, candi, covars, use.coords = FALSE, strata.type = "area", 
            x.max, x.min, y.max, y.min, iterations,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            boundary, progress = TRUE, verbose = TRUE, greedy = FALSE,
            weights, nadir, utopia) {
    
    if (!is.data.frame(covars)) covars <- as.data.frame(covars)
    
    # Check spsann arguments ###################################################
    eval(.check_spsann_arguments())
    ############################################################################
    
    check <- .optimCORRcheck(candi = candi, covars = covars, 
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Set plotting options ####################################################
    eval(.plotting_options())
    ############################################################################
    
    # Prepare sample points
    #n_candi <- nrow(candi)
    #points <- .spsannPoints(points = points, candi = candi, n.candi = n_candi)
    #n_pts <- nrow(points)
    #conf0 <- points
    #old_conf <- conf0
    # Prepare points and candi##################################################
    prepare_points <- 
      function (...) {parse(text = readLines("tools/prepare-points.R"))}
    eval(prepare_points())
    ############################################################################
    
    # Prepare for jittering ####################################################
    prepare_jittering <- 
      function (...) {parse(text = readLines("tools/prepare-jittering.R"))}
    eval(prepare_jittering())
    ############################################################################
    
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
      scm <- cor(sm, use = "complete.obs")
      energy0 <- sum(abs(pcm - scm))
    } else { # Factor covariates
      if (covars.type == "factor") {
        pcm <- pedometrics::cramer(covars)
        scm <- pedometrics::cramer(sm)
        energy0 <- sum(abs(pcm - scm))
      }
    }

    # Other settings for the simulated annealing algorithm
    MOOP <- FALSE
    old_scm      <- scm
    new_scm      <- scm
    best_scm     <- scm
    old_sm       <- sm
    new_sm       <- sm
    best_sm      <- sm
    count        <- 0
    old_energy   <- energy0
    best_energy  <- Inf
    energies     <- vector()
    accept_probs <- vector()
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()

    # Begin the main loop
    for (k in 1:iterations) {

      # Jitter one of the points and update x.max and y.max; which point (wp)?
      wp <- sample(c(1:n_pts), 1)
      new_conf <- spJitterFinite(old_conf, candi, x.max, x.min, y.max, 
                                 y.min, wp)
      x.max <- x_max0 - (k / iterations) * (x_max0 - x.min)
      y.max <- y_max0 - (k / iterations) * (y_max0 - y.min)

      # Update sample and correlation matrices, and energy state
      if (covars.type == "numeric") { # Numeric covariates
        new_row <- covars[new_conf[wp, 1], ]
        new_sm[wp, ] <- new_row
        new_scm <- cor(new_sm, use = "complete.obs")
        new_energy <- sum(abs(pcm - new_scm))

      } else { # Factor covariates
        if (covars.type == "factor") {
          new_row <- covars[new_conf[wp, 1], ]
          new_sm[wp, ] <- new_row
          new_scm <- pedometrics::cramer(new_sm)
          new_energy <- sum(abs(pcm - new_scm))
        }
      }
      
      # Evaluate the new system configuration
      if (greedy) {
        random_prob <- 1
      } else {
        random_prob <- runif(1)
      }
      actual_prob <- acceptance[[1]] * exp(-k / acceptance[[2]])
      accept_probs[k] <- actual_prob
      if (new_energy <= old_energy) {
        old_conf   <- new_conf
        old_energy <- new_energy
        count      <- 0
        old_sm     <- new_sm
        old_scm    <- new_scm
      } else {
        if (new_energy > old_energy & random_prob <= actual_prob) {
          old_conf   <- new_conf
          old_energy <- new_energy
          count      <- count + 1
          old_sm     <- new_sm
          old_scm    <- new_scm
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... p = ",
                random_prob, "\n")
          }
        } else {
          new_energy <- old_energy
          new_conf   <- old_conf
          count      <- count + 1
          new_sm     <- old_sm
          new_scm    <- old_scm
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... stops at",
                stopping[[1]], "\n")
          }
        }
      }
      # Best energy state
      energies[k] <- new_energy
      if (new_energy < best_energy / 1.0000001) {
        best_k          <- k
        best_conf       <- new_conf
        best_energy     <- new_energy
        best_old_energy <- old_energy
        old_conf        <- old_conf
        best_sm         <- new_sm
        best_old_sm     <- old_sm
        best_scm        <- new_scm
        best_old_scm    <- old_scm
      }
      
      # Plotting
      if (plotit && pedometrics::is.numint(k / 10)) {
        .spSANNplot(energy0 = energy0, energies = energies, k = k, 
                    acceptance = acceptance, accept_probs = accept_probs,
                    boundary = boundary, new_conf = new_conf[, 2:3], 
                    conf0 = conf0[, 2:3], y_max0 = y_max0, y.max = y.max, 
                    x_max0 = x_max0, x.max = x.max, best.energy = best_energy,
                    best.k = best_k, MOOP = MOOP)
      }
      
      # Freezing parameters
      if (count == stopping[[1]]) {
        if (new_energy > best_energy * 1.000001) {
          old_conf   <- old_conf
          new_conf   <- best_conf
          old_energy <- best_old_energy
          new_energy <- best_energy
          count      <- 0
          new_sm     <- best_sm
          new_scm    <- best_scm
          old_sm     <- best_old_sm
          old_scm    <- best_old_scm
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
    res <- .spSANNout(new_conf = new_conf, energy0 = energy0, 
                      energies = energies, time0 = time0, MOOP = MOOP)
    return (res)
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
.optimCORRcheck <-
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
    st <- match(strata.type, c("area", "range"))
    if (is.na(st)) {
      res <- paste("'strata.type = ", strata.type, "' is not supported",
                   sep = "")
      return (res)
    }
  }
# FUNCTION - CALCULATE ENERGY STATE ############################################
#' @rdname optimCORR
#' @export
objCORR <-
  function (points, candi, covars, use.coords = FALSE, strata.type = "area") {
    
    if (!is.data.frame(covars)) covars <- as.data.frame(covars)
    
    # Check arguments
    check <- .optimCORRcheck(candi = candi, covars = covars, 
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare sample points
    #n_candi <- nrow(candi)
    #points <- .spsannPoints(points = points, candi = candi, n.candi = n_candi)
    #n_pts <- nrow(points)
    # Prepare points and candi##################################################
    prepare_points <- 
      function (...) {parse(text = readLines("tools/prepare-points.R"))}
    eval(prepare_points())
    ############################################################################
    
    # Prepare covariates (covars) and create the starting sample matrix (sm)
    #     if (use.coords) {
    #       covars <- .useCoords(covars.type = covars.type, candi = candi, 
    #                            n.pts = n_pts, strata.type = strata.type)
    #     }
    #     sm <- covars[points[, 1], ]
    covars.type <- ifelse(is.any.factor(covars), "factor", "numeric")
    covars <- .covarsACDC(covars = covars, covars.type = covars.type, 
                          use.coords = use.coords, candi = candi, n.pts = n_pts,
                          strata.type = strata.type)
    sm <- covars[points[, 1], ]
    
    # Calculate the energy state
    if (covars.type == "numeric") { # Numeric covariates
      pcm <- cor(covars, use = "complete.obs")
      scm <- cor(sm, use = "complete.obs")
      energy <- sum(abs(pcm - scm))
    } else { # Factor covariates
      if (covars.type == "factor") {
        pcm <- pedometrics::cramer(covars)
        scm <- pedometrics::cramer(sm)
        energy <- sum(abs(pcm - scm))
      }
    }
    return (energy)
  }
