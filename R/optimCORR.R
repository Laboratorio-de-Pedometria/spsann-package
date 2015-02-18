#' Optimization of sample patterns for trend estimation
#'
#' Optimize a sample pattern for trend estimation. The criterion is defined so 
#' that the sample reproduces the association/correlation between the covariates
#' (\bold{CORR}).
#'
#' @template spJitter_doc
#' @template spSANN_doc
#'
#' @param covars Data frame or matrix with the covariates in the columns.
#'
#' @param covars.type Character value. The type of covariates that is
#' being used. Available options are \code{"numeric"} and \code{"factor"}.
#' Defaults to \code{covars.type = "numeric"}.
#'
#' @param use.coords Logical value. Should the coordinates be used as 
#' covariates? Defaults to \code{use.coords = FALSE}.
#'
#' @param strata.type Character value. The type of strata to be used to 
#' categorize the coordinates when they are used with covariates of type factor.
#' Available options are \code{"area"} for equal area and \code{"range"} for
#' equal range. Defaults to \code{strata.type = "area"}.
#' 
#' @details
#' This method was derived from the conditioned Latin Hypercube of Minasny and
#' McBratney (2006). Visit the package manual to see the corrections that we
#' have made in that method.
#'
#' @return
#' \code{optimCORR} returns a matrix: the optimized sample pattern with
#' the evolution of the energy state during the optimization as an attribute.
#' 
#' \code{objCORR} returns a numeric value: the energy state of the point 
#' pattern.
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
#' @importFrom pedometrics cramer
#' @importFrom pedometrics is.numint
#' @importFrom pedometrics cont2cat
#' @importFrom SpatialTools dist2
#' @export
#' @examples
#' require(pedometrics)
#' require(sp)
#' require(rgeos)
#' data(meuse.grid)
#' candi              <- meuse.grid[, 1:2]
#' coordinates(candi) <- ~ x + y
#' gridded(candi)     <- TRUE
#' boundary           <- as(candi, "SpatialPolygons")
#' boundary           <- gUnionCascaded(boundary)
#' candi              <- coordinates(candi)
#' candi              <- matrix(cbind(1:dim(candi)[1], candi), ncol = 3)
#' covars             <- meuse.grid[, 5]
#' x.max              <- diff(bbox(boundary)[1, ])
#' y.max              <- diff(bbox(boundary)[2, ])
#' y.min              <- 40
#' x.min              <- 40
#' set.seed(2001)
#' res <- optimCORR(points = 100, candi = candi, covars = covars, 
#'                  use.coords = TRUE, covars.type = "numeric", x.max = x.max, 
#'                  x.min = x.min, y.max = y.max, y.min = y.min,
#'                  boundary = boundary, iterations = 500)
#' tail(attr(res, "energy"))
#' objCORR(points = res, candi = candi, covars = covars, 
#'         covars.type = "numeric", use.coords = TRUE)
# MAIN FUNCTION ################################################################
optimCORR <-
  function (points, candi, covars, covars.type = "numeric", use.coords = FALSE, 
            strata.type = "area", x.max, x.min, y.max, y.min, iterations,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            boundary, progress = TRUE, verbose = TRUE) {
    
    if (!is.data.frame(covars)) covars <- as.data.frame(covars)
    
    # Check arguments
    check <- .spSANNcheck(points, candi, x.max, x.min, y.max, y.min,
                          iterations, acceptance, stopping, plotit, boundary,
                          progress, verbose)
    if (!is.null(check)) stop (check, call. = FALSE)
    check <- .optimCORRcheck(candi = candi, covars = covars, 
                             covars.type = covars.type, use.coords = use.coords,
                             strata.type = strata.type)
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
    if (use.coords) {
      if (covars.type == "factor") {
        coords <- data.frame(candi[, 2:3])
        breaks <- .coordStrata(n_pts, coords, strata.type)
        coords <- pedometrics::cont2cat(coords, breaks)
        covars <- data.frame(covars, coords)
      } else {
        covars <- data.frame(covars, candi[, 2:3])
      }
    }
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
      random_prob <- runif(1)
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
      if (plotit && any(round(seq(1, iterations, 10)) == k)) {
        .spSANNplot(energy0, energies, k, acceptance,
                    accept_probs, boundary, new_conf[, 2:3],
                    conf0[, 2:3], y_max0, y.max, x_max0, x.max)
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
    res <- .spSANNout(new_conf, energy0, energies, time0)
    return (res)
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
.optimCORRcheck <-
  function (candi, covars, covars.type, use.coords, strata.type) {
    
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
    if (missing(covars.type)) {
      res <- paste("'covars.type' is missing")
      return (res)
    } else {
      ct <- pmatch(covars.type, c("numeric", "factor"))
      if (is.na(ct)) {
        res <- paste("'covars.type = ", covars.type, "' is not supported", 
                     sep = "")
        return (res)
      }
    }
    
    # strata.type
    st <- match(strata.type, c("area", "range"))
    if (is.na(st)) {
      res <- paste("'strata.type = ", strata.type, "' is not supported",
                   sep = "")
      return (res)
    }
  }
# INTERNAL FUNCTION - BREAKS FOR COORDINATES ###################################
.coordStrata <-
  function (n.pts, coords, strata.type) {
    # equal area strata
    if (strata.type == "area") {
      n_cov <- 2
      probs <- seq(0, 1, length.out = n.pts + 1)
      breaks <- lapply(coords, quantile, probs, na.rm = TRUE, type = 3)
      breaks <- lapply(breaks, unique)
    } else {
      # equal range strata
      if (strata.type == "range") {
        n_cov <- 2
        breaks <- lapply(1:n_cov, function(i)
          seq(min(covars[, i]), max(covars[, i]), length.out = n.pts + 1))
        d <- lapply(1:n_cov, function(i)
          SpatialTools::dist2(matrix(breaks[[i]]), matrix(covars[, i])))
        d <- lapply(1:n_cov, function(i) apply(d[[i]], 1, which.min))
        breaks <- lapply(1:n_cov, function(i) breaks[[i]] <- covars[d[[i]], i])
        breaks <- lapply(breaks, unique)
      }
    }
    return (breaks)
  }
# FUNCTION - CALCULATE ENERGY STATE ############################################
#' @rdname optimCORR
#' @export
objCORR <-
  function (points, candi, covars, covars.type = "numeric", use.coords = FALSE, 
            strata.type = "area") {
    
    if (!is.data.frame(covars)) covars <- as.data.frame(covars)
    
    # Check arguments
    check <- .optimCORRcheck(candi = candi, covars.type = covars.type,
                             covars = covars, use.coords = use.coords,
                             strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare sample points
    n_candi <- nrow(candi)
    points <- .spsannPoints(points = points, candi = candi, n.candi = n_candi)
    n_pts <- nrow(points)
    
    # Prepare covariates (covars) and create the starting sample matrix (sm)
    if (use.coords) {
      if (covars.type == "factor") {
        coords <- data.frame(candi[, 2:3])
        breaks <- .coordStrata(n_pts, coords, strata.type)
        coords <- pedometrics::cont2cat(coords, breaks)
        covars <- data.frame(covars, coords)
      } else {
        covars <- data.frame(covars, candi[, 2:3])
      }
    }
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
