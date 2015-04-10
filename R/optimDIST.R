#' Optimization of sample patterns for trend estimation
#'
#' Optimize a sample pattern for trend estimation. A criterion is defined so 
#' that the sample reproduces the marginal distribution of the covariates
#' (\bold{DIST}).
#'
#' @template spJitter_doc
#' @template spSANN_doc
#'
#' @param covars Data frame or matrix with the covariates in the columns.
#'
#' @param use.coords Logical value. Should the coordinates be used as 
#' covariates? Defaults to \code{use.coords = FALSE}.
#'
#' @param strata.type Character value. The type of strata to be used with
#' numeric covariates. Available options are \code{"area"} for equal area and
#' \code{"range"} for equal range. Defaults to \code{strata.type = "area"}. See
#' \sQuote{Details} for more information.
#' 
#' @param greedy Logical value. Should the optimization be done using a greedy
#' algorithm, that is, without accepting worse system configurations? Defaults
#' to \code{greedy = FALSE}.
#' 
#' @details
#' This method derives from the method known as the conditioned Latin Hypercube
#' originally proposed by Minasny and McBratney (2006). Visit the package manual
#' to see the improvements that we have made in that method.
#'
#' @return
#' \code{optimDIST} returns a matrix: the optimized sample pattern with
#' the evolution of the energy state during the optimization as an attribute.
#'  
#' \code{objDIST} returns a numeric value: the energy state of the sample
#' pattern - the objective function value.
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
#' algorithm incorporating operational constraints. \emph{5th Global Workshop on
#' Digital Soil Mapping}. Sydney, p. 227-231, 2012.
#'
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[clhs]{clhs}}
#' @keywords spatial optimize
#' @concept simulated annealing
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
#' candi <- matrix(cbind(1:dim(candi)[1], candi), ncol = 3)
#' covars <- meuse.grid[, 5]
#' x.max <- diff(bbox(boundary)[1, ])
#' y.max <- diff(bbox(boundary)[2, ])
#' set.seed(2001)
#' res <- optimDIST(points = 100, candi = candi, covars = covars, 
#'                  use.coords = TRUE, x.max = x.max, x.min = 40, y.max = y.max, 
#'                  y.min = 40, boundary = boundary, iterations = 1000)
#' tail(attr(res, "energy"), 1) # 0.9897776
#' objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)
# MAIN FUNCTION ################################################################
optimDIST <-
  function (points, candi, covars, strata.type = "area", use.coords = FALSE, 
            x.max, x.min, y.max, y.min, iterations,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            boundary, progress = TRUE, verbose = TRUE, greedy = FALSE) {
    
    if (!is.data.frame(covars)) covars <- as.data.frame(covars)
    
    # Check arguments
    check <- .spSANNcheck(points, candi, x.max, x.min, y.max, y.min,
                          iterations, acceptance, stopping, plotit, boundary,
                          progress, verbose)
    if (!is.null(check)) stop (check, call. = FALSE)
    check <- .optimDISTcheck(candi = candi, covars = covars, 
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    if (plotit) {
      par0 <- par()
      on.exit(suppressWarnings(par(par0)))
    }
    
    # Prepare sample points
    n_candi <- nrow(candi)
    points <- .spsannPoints(points = points, candi = candi, n.candi = n_candi)
    n_pts <- nrow(points)
    conf0 <- points
    old_conf <- conf0
    
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
                    best.k = best_k)
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
    res <- .spSANNout(new_conf, energy0, energies, time0)
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
    n_candi <- nrow(candi)
    points <- .spsannPoints(points = points, candi = candi, n.candi = n_candi)
    n_pts <- nrow(points)
    
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
