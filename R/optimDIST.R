#' Optimization of sample configurations for spatial trend estimation
#'
#' Optimize a sample configuration for trend estimation. A criterion is defined 
#' so that the sample reproduces the marginal distribution of the covariates
#' (\bold{DIST}).
#'
#' @template spJitter_doc
#' @template spSANN_doc
#' @template ACDC_doc
#' 
#' @details
#' Reproducing the marginal distribution of the numeric covariates depends upon
#' the definition of sampling strata. These sampling strata are also used to 
#' categorize any numeric covariate when they are passed together with 
#' categorical covariates (see more details at \code{optimDist}).
#' Two types of sampling strata can be used. \emph{Equal-area} sampling strata 
#' are defined using the sample quantiles estimated with \code{quantile()} using
#' a discontinuous function (\code{type = 3}). This is to avoid creating 
#' breakpoints that do not occur in the population of existing covariate values.
#' 
#' The function \code{quantile()} commonly produces repeated break points. A
#' break point will always be repeated if that value has a relatively
#' high frequency in the population of covariate values. The number of repeated
#' break points increases with the number of sampling strata. Only unique
#' break points are used to create sampling strata.
#' 
#' \emph{Equal-range} sampling strata are defined breaking the range of 
#' covariate values into pieces of equal size. This method usually creates
#' breakpoints that do not occur in the population of existing covariate values.
#' Such breakpoints are replaced by the nearest existing covariate value 
#' identified using Euclidean distances.
#' 
#' Both stratification methods can produce sampling strata that cover a range of
#' values that do not exist in the population of covariate value. Any empty 
#' sampling strata is merged with the closest non-empty sampling strata. These
#' are identified using Euclidean distances.
#' 
#' The approaches used to define the sampling strata result in each numeric 
#' covariate having a different number of sampling strata, some of them with 
#' different area/size. Because the goal is to have a sample that reproduces the
#' marginal distribution of the covariate, each sampling strata will have a
#' different number of sample points. The wanted distribution of the number of 
#' sample points per strata is estimated empirically computing the proportion of
#' points of the population of existing covariate values that fall in each
#' sampling strata.
#' 
#' @return
#' \code{optimDIST} returns a matrix: the optimized sample configuration with
#' the evolution of the energy state during the optimization as an attribute.
#'  
#' \code{objDIST} returns a numeric value: the energy state of the sample
#' configuration - the objective function value.
#'
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[clhs]{clhs}}
#' @importFrom pedometrics is.numint
#' @importFrom pedometrics cont2cat
#' @importFrom pedometrics is.all.factor
#' @importFrom pedometrics is.any.factor
#' @importFrom SpatialTools dist2
#' @import Rcpp
#' @export
#' @examples
#' require(pedometrics)
#' require(sp)
#' require(rgeos)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' coordinates(candi) <- ~ x + y
#' gridded(candi) <- TRUE
#' boundary <- as(candi, "SpatialPolygons")
#' boundary <- gUnionCascaded(boundary)
#' candi <- coordinates(candi)
#' covars <- meuse.grid[, 5]
#' x.max <- diff(bbox(boundary)[1, ])
#' y.max <- diff(bbox(boundary)[2, ])
#' set.seed(2001)
#' res <- optimDIST(points = 100, candi = candi, covars = covars, 
#'                  use.coords = TRUE, x.max = x.max, x.min = 40, y.max = y.max, 
#'                  y.min = 40, boundary = boundary, iterations = 100)
#' tail(attr(res, "energy"), 1) # 1.656926
#' objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)
# MAIN FUNCTION ################################################################
optimDIST <-
  function (points, candi, covars, strata.type = "area", use.coords = FALSE, 
            x.max, x.min, y.max, y.min, iterations,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            boundary, progress = TRUE, verbose = TRUE, greedy = FALSE,
            weights, nadir, utopia) {
    
    if (!is.data.frame(covars)) covars <- as.data.frame(covars)
    
    # Check spsann arguments ###################################################
    check_spsann_arguments <- 
      function (...) {parse(text = readLines("tools/check-spsann-arguments.R"))}
    eval(check_spsann_arguments())
    ############################################################################
    check <- .optimDISTcheck(candi = candi, covars = covars, 
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    if (plotit) {
      par0 <- par()
      on.exit(suppressWarnings(par(par0)))
    }
    
    # Prepare sample points
    #n_candi <- nrow(candi)
    #points <- .spsannPoints(points = points, candi = candi, n.candi = n_candi)
    #n_pts <- nrow(points)
    #conf0 <- points
    #old_conf <- conf0
    # Prepare points and candi #################################################
    prepare_points <- 
      function (...) {parse(text = readLines("tools/prepare-points.R"))}
    eval(prepare_points())
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
      strata <- .numStrata(n.pts = n_pts, covars = covars, 
                           strata.type = strata.type)
      energy0 <- .objDISTnum(sm = sm, n.pts = n_pts, n.cov = n_cov, 
                             strata = strata)
    } else { # Factor covariates
      if (covars.type == "factor") {
        #pop_prop <- lapply(covars, function(x) table(x) / nrow(covars) * 100)
        pop_prop <- lapply(covars, function(x) table(x) / nrow(covars))
        energy0 <- .objDISTfac(sm = sm, pop.prop = pop_prop, n.pts = n_pts,
                                n.cov = n_cov)
      }
    }

    # Other settings for the simulated annealing algorithm
    MOOP <- FALSE
    old_sm       <- sm
    new_sm       <- sm
    best_sm      <- sm
    count        <- 0
    old_energy   <- energy0
    best_energy  <- Inf
    energies     <- vector()
    accept_probs <- vector()
    x_max0       <- x.max
    y_max0       <- y.max
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
        new_energy <- .objDISTnum(sm = new_sm, n.cov = n_cov, strata = strata, 
                                  n.pts = n_pts)
      } else { # Factor covariates
        if (covars.type == "factor") {
          new_row <- covars[new_conf[wp, 1], ]
          new_sm[wp, ] <- new_row
          new_energy <- .objDISTfac(sm = new_sm, pop.prop = pop_prop, 
                                    n.pts = n_pts, n.cov = n_cov)
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
      } else {
        if (new_energy > old_energy & random_prob <= actual_prob) {
          old_conf   <- new_conf
          old_energy <- new_energy
          count      <- count + 1
          old_sm     <- new_sm
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... p = ",
                random_prob, "\n")
          }
        } else {
          new_energy <- old_energy
          new_conf   <- old_conf
          count      <- count + 1
          new_sm     <- old_sm
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
      }
      # Plotting
      #if (plotit && any(round(seq(1, iterations, 10)) == k)) {
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
          old_sm     <- best_old_sm
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
    res <- .spSANNout(new_conf, energy0, energies, time0, MOOP = MOOP)
    return (res)
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
.optimDISTcheck <-
  function (candi, covars, 
            #covars.type, 
            use.coords, strata.type) {
    
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
    
    # covars.type
#     if (missing(covars.type)) {
#       res <- paste("'covars.type' is missing")
#       return (res)
#     } else {
#       ct <- pmatch(covars.type, c("numeric", "factor"))
#       if (is.na(ct)) {
#         res <- paste("'covars.type = ", covars.type, "' is not supported", 
#                      sep = "")
#         return (res)
#       }
#     }
    
    # strata.type
    st <- match(strata.type, c("area", "range"))
    if (is.na(st)) {
      res <- paste("'strata.type = ", strata.type, "' is not supported",
                   sep = "")
      return (res)
    }
  }
# INTERNAL FUNCTION - CRITERION FOR FACTOR COVARIATES ##########################
.objDISTfac <-
  function (sm, pop.prop, n.pts, n.cov) {    
    #samp_prop <- lapply(sm, function(x) table(x) / n.pts * 100)
    samp_prop <- lapply(sm, function(x) table(x) / n.pts)
    samp_prop <- sapply(1:n.cov, function (i)
      sum(abs(samp_prop[[i]] - pop.prop[[i]])))
    energy <- sum(samp_prop)
    return (energy)
  }
# INTERNAL FUNCTION - CRITERION FOR NUMERIC COVARIATES #########################
.objDISTnum <-
  function (sm, n.pts, n.cov, strata) {
    counts <- lapply(1:n.cov, function (i)
      hist(sm[, i], strata[[1]][[i]], plot = FALSE)$counts)
    #counts <- lapply(1:n.cov, function(i) counts[[i]] / n.pts * 100)
    counts <- lapply(1:n.cov, function(i) counts[[i]] / n.pts)
    counts <- sapply(1:n.cov, function (i) 
      sum(abs(counts[[i]] - strata[[2]][[i]])))
    energy <- sum(counts)
    return (energy)
  }
# CALCULATE OBJECTIVE FUNCTION VALUE ###########################################
#' @rdname optimDIST
#' @export
objDIST <-
  function (points, candi, covars, strata.type = "area", use.coords = FALSE) {
    
    if (!is.data.frame(covars)) covars <- as.data.frame(covars)
    
    # Check arguments
    check <- .optimDISTcheck(candi = candi, covars = covars,
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare sample points
    #n_candi <- nrow(candi)
    #points <- .spsannPoints(points = points, candi = candi, n.candi = n_candi)
    #n_pts <- nrow(points)
    # Prepare points and candi #################################################
    prepare_points <- 
      function (...) {parse(text = readLines("tools/prepare-points.R"))}
    eval(prepare_points())
    ############################################################################
    
    # Prepare covariates (covars) and create the starting sample matrix (sm)
    covars.type <- ifelse(is.any.factor(covars), "factor", "numeric")
    covars <- .covarsACDC(covars = covars, covars.type = covars.type, 
                          use.coords = use.coords, candi = candi, n.pts = n_pts,
                          strata.type = strata.type)
    n_cov <- ncol(covars)
    sm <- covars[points[, 1], ]
    
    # Calculate the energy state
    if (covars.type == "numeric") { # Numeric covariates
      strata <- .numStrata(n.pts = n_pts, covars = covars, 
                           strata.type = strata.type)
      energy <- .objDISTnum(sm = sm, n.pts = n_pts, n.cov = n_cov, 
                            strata = strata)
    } else { # Factor covariates
      if (covars.type == "factor") {
        #pop_prop <- lapply(covars, function(x) table(x) / nrow(covars) * 100)
        pop_prop <- lapply(covars, function(x) table(x) / nrow(covars))
        energy <- .objDISTfac(sm = sm, pop.prop = pop_prop, n.pts = n_pts,
                              n.cov = n_cov)
      }
    }
    return (energy)
  }
