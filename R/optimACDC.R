#' Optimization of sample patterns for trend estimation
#'
#' Optimize a sample pattern for trend estimation. The criterion is defined so 
#' that the sample reproduces the association/correlation between the
#' covariates, as well as their marginal distribution (\bold{ACDC}).
#'
#' @template spJitter_doc
#' @template spSANN_doc
#'
#' @param covars Data frame or matrix with the covariates in the columns.
#'
#' @param weights List with two named sub-arguments. The weights assigned to 
#' correlation/association measure and the sampling strata/classes. They
#' must be named and sum to unity. Defaults to 
#' \code{weights = list(correl = 0.5, strata = 0.5)}.
#'
#' @param use.coords Logical value. Should the coordinates be used as 
#' covariates? Defaults to \code{use.coords = FALSE}.
#'
#' @param strata.type Character value. The type of strata to be used with
#' numeric covariates. Available options are \code{"area"} for equal area and
#' \code{"range"} for equal range. Defaults to \code{strata.type = "area"}. See
#' \sQuote{Details} for more information.
#'
#' @param nadir List with four named sub-arguments: \code{sim} -- the number of
#' random realizations to estimate the nadir point; \code{save.sim} -- logical 
#' for saving the simulated values and returning them as an attribute of the 
#' optimized sample pattern; \code{user} -- a user-defined value;
#' \code{abs} -- logical for calculating the nadir point internally. Only
#' simulations are implemented in the current version. Defaults to 
#' \code{nadir = list(sim = 1000, save.sim = TRUE, user = NULL, abs = NULL)}.
#' 
#' @param utopia List with two named sub-arguments: \code{user} -- a list of
#' two user-defined values (\code{correl} and \code{strata}), and \code{abs} --
#' a logical value for calculating the utopia point internally. Defaults to 
#' \code{utopia = list(user = list(correl = NULL, strata = NULL), abs = NULL)}.
#' 
#' @param scale List with two named sub-arguments: \code{type} -- the type of
#' scaling that should be used, with available options \code{"upper} and
#' \code{"upper-lower}; and \code{max} -- the maximum value after scaling.
#' Defaults to \code{scale = list(type = "upper-lower", max = 100)}.
#' 
#' @details
#' This method is also known as the conditioned Latin Hypercube of Minasny and
#' McBratney (2006). Visit the package manual to see the corrections that we
#' have made in that method.
#'
#' @return
#' \code{optimACDC} returns a matrix: the optimized sample pattern with
#' the evolution of the energy state during the optimization as an attribute.
#' 
#' \code{objACDC} returns a numeric value: the energy state of the point 
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
#' @importFrom pedometrics is.all.factor
#' @importFrom pedometrics is.any.factor
#' @importFrom SpatialTools dist2
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
#' y.min <- 40
#' x.min <- 40
#' y.max <- diff(bbox(boundary)[2, ])
#' nadir <- list(sim = 10, save.sim = TRUE, user = NULL, abs = NULL)
#' utopia <- list(user = list(correl = 0, strata = 0), abs = NULL)
#' scale <- list(type = "upper-lower", max = 100)
#' weights <- list(strata = 0.5, correl = 0.5)
#' set.seed(2001)
#' res <- optimACDC(points = 100, candi = candi, covars = covars,
#'                  use.coords = TRUE, covars.type = "numeric", 
#'                  weights = weights, x.max = x.max, x.min = x.min, 
#'                  y.max = y.max, y.min = y.min, boundary = boundary, 
#'                  nadir = nadir, iterations = 500, utopia = utopia, 
#'                  scale = scale)
#' tail(attr(res, "energy"))
#' objACDC(points = res, candi = candi, covars = covars, use.coords = TRUE, 
#'         covars.type = "numeric", weights = weights, nadir = nadir,
#'         utopia = utopia, scale = scale)
# MAIN FUNCTION ################################################################
optimACDC <-
  function (points, candi, covars, strata.type = "area", 
            weights = list(correl = 0.5, strata = 0.5), use.coords = FALSE,
            nadir = list(sim = 1000, save.sim = TRUE, user = NULL, abs = NULL),
            utopia = list(user = list(correl = NULL, strata = NULL), abs = NULL), 
            scale = list(type = "upper-lower", max = 100),
            x.max, x.min, y.max, y.min, iterations,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            boundary, progress = TRUE, verbose = TRUE) {
    
    if (!is.data.frame(covars)) covars <- as.data.frame(covars)
    
    # Check arguments
    check <- .spSANNcheck(points, candi, x.max, x.min, y.max, y.min,
                          iterations, acceptance, stopping, plotit, boundary,
                          progress, verbose)
    if (!is.null(check)) stop (check, call. = FALSE)
    check <- .optimACDCcheck(candi = candi, covars = covars, nadir = nadir,
                             weights = weights, use.coords = use.coords, 
                             strata.type = strata.type, utopia = utopia, 
                             scale = scale)
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
      pcm <- cor(covars, use = "complete.obs")
      strata <- .numStrata(n.pts = n_pts, covars = covars, 
                           strata.type = strata.type)
      nadir <- .numNadir(n.pts = n_pts, n.cov = n_cov, n.candi = n_candi, 
                         pcm = pcm, nadir = nadir, candi = candi, 
                         covars = covars, strata = strata, scale = scale)
      utopia <- .numUtopia(utopia = utopia)
      scm <- cor(sm, use = "complete.obs")
      energy0 <- .objNum(sm = sm, n.cov = n_cov, strata = strata, pcm = pcm, 
                         scm = scm, nadir = nadir, weights = weights, 
                         n.pts = n_pts, utopia = utopia, scale = scale)

    } else { # Factor covariates
      if (covars.type == "factor") {
        pcm <- pedometrics::cramer(covars)
        # ASR: We multiply the proportions by 100 for numerical stability
        pop_prop <- lapply(covars, function(x) table(x) / nrow(covars))
        #pop_prop <- lapply(covars, function(x) table(x) / nrow(covars) * 100)
        nadir <- .facNadir(nadir = nadir, candi = candi, n.candi = n_candi,
                           n.pts = n_pts, n.cov = n_cov, covars = covars, 
                           pop.prop = pop_prop, pcm = pcm, scale = scale)
        utopia <- .facUtopia(utopia = utopia)
        scm <- pedometrics::cramer(sm)
        energy0 <- .objFac(sm = sm, pop.prop = pop_prop, nadir = nadir, 
                           weights = weights, pcm = pcm, scm = scm,
                           n.pts = n_pts, n.cov = n_cov, utopia = utopia,
                           scale = scale)
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
        new_energy <- .objNum(sm = new_sm, n.cov = n_cov, strata = strata, 
                              pcm = pcm, scm = new_scm, nadir = nadir,
                              weights = weights, n.pts = n_pts, utopia = utopia,
                              scale = scale)

      } else { # Factor covariates
        if (covars.type == "factor") {
          new_row <- covars[new_conf[wp, 1], ]
          new_sm[wp, ] <- new_row
          new_scm <- pedometrics::cramer(new_sm)
          new_energy <- .objFac(sm = new_sm, pop.prop = pop_prop, scm = new_scm, 
                                nadir = nadir, weights = weights, pcm = pcm, 
                                n.pts = n_pts, n.cov = n_cov, utopia = utopia,
                                scale = scale)
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
.optimACDCcheck <-
  function (candi, covars, weights, use.coords, strata.type, nadir, utopia,
            scale) {
    
    # utopia
    aa <- !is.list(utopia)
    bb <- length(utopia) != 2
    cc <- is.null(names(utopia))
    dd <- !all(c(names(utopia) == c("user", "abs")) == TRUE)
    if (aa || bb || cc || dd) {
      res <- paste("'utopia' must be a list with two named sub-arguments:",
                   "'user' and 'abs'", sep = "")
      return (res)
    }
    
    # scale
    aa <- !is.list(scale)
    bb <- length(scale) != 2
    cc <- is.null(names(scale))
    dd <- !all(c(names(scale) == c("type", "max")) == TRUE)
    if (aa || bb || cc || dd) {
      res <- paste("'scale' must be a list with two named sub-arguments:",
                   "'type' and 'max'", sep = "")
      return (res)
    }
    if (!match(scale$type, c("upper", "upper-lower"))) {
      res <- paste("'scale$type' must be 'upper' or 'upper-lower'", sep = "")
      return (res)
    }
    
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
    
    # weights
    aa <- !is.list(weights)
    bb <- length(weights) != 2
    cc <- is.null(names(weights))
    #dd <- !all(c(names(weights) == c("correl", "strata")) == TRUE)
    #if (aa || bb || cc || dd) {
    if (aa || bb || cc || dd) {
      res <- paste("'weights' must be a list with two named sub-arguments:",
                   "'strata' and 'correl'", sep = "")
      return (res)
    }
    if (sum(unlist(weights)) != 1) {
      res <- paste("the 'weights' must sum to 1")
      return (res)
    }
    
    # strata.type
    st <- match(strata.type, c("area", "range"))
    if (is.na(st)) {
      res <- paste("'strata.type = ", strata.type, "' is not supported",
                   sep = "")
      return (res)
    }
    
    # nadir
    if (!is.list(nadir) || length(nadir) != 4) {
      res <- paste("'nadir' must be a list with four sub-arguments")
      return (res)
    }
    n <- !sapply(nadir, is.null)
    if (n[[1]] == TRUE) {
      if (n[[2]] == FALSE) {
        res <- paste("you must inform if the simulations should be saved")
        return (res)
      }
      if (n[[3]] == TRUE || n[[4]] == TRUE) {
        res <- paste("you must choose a single nadir option")
        return (res)
      }
    } else {
      if (n[[3]] == TRUE) {
        res <- paste("sorry but you cannot set the nadir point")
        return (res)
      }
      if (n[[4]] == TRUE) {
        res <- paste("sorry but the nadir point cannot be calculated")
        return (res)
      }
    }
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
      
      breaks <- lapply(breaks, unique)
      count <- lapply(1:n_cov, function (i)
        hist(covars[, i], breaks[[i]], plot = FALSE)$counts)
      
      # ASR: We use the proportion of points per sampling strata
      #count <- lapply(1:n_cov, function(i) 
      #  count[[i]] / sum(count[[i]]) * n.pts)
      
      #count <- lapply(1:n_cov, function(i) count[[i]] / sum(count[[i]]) * 100)
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
        
        # ASR: We use the proportion of points per sampling strata
        #count <- lapply(1:n_cov, function(i)
        #  count[[i]] / sum(count[[i]]) * n.pts)
        count <- lapply(1:n_cov, function(i)
          #count[[i]] / sum(count[[i]]) * 100)
          count[[i]] / sum(count[[i]]))
        
        
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
  function (n.pts, n.cov, n.candi, pcm, nadir, candi, covars, strata, scale) {

    if (!is.null(nadir[[1]]) && !is.null(nadir[[2]])) { # Simulate the nadir
      m <- paste("simulating ", nadir[[1]], " nadir values...", sep = "")
      message(m)
      
      # SET VARIABLES
      strata_nadir <- vector()
      correl_nadir <- vector()
      
      # BEGIN THE SIMULATION
      for (i in 1:nadir[[1]]) {
        pts <- sample(1:n.candi, n.pts)
        sm <- covars[pts, ]
        scm <- cor(sm, use = "complete.obs")
        counts <- lapply(1:n.cov, function (i)
          hist(sm[, i], strata[[1]][[i]], plot = FALSE)$counts)
        
        # ASR: We use the proportion of points per sampling strata
        #counts <- lapply(1:n.cov, function(i) counts[[i]] / n.pts * 100)
        counts <- lapply(1:n.cov, function(i) counts[[i]] / n.pts)
        counts <- sapply(1:n.cov, function (i)
          sum(abs(counts[[i]] - strata[[2]][[i]])))
        strata_nadir[i] <- sum(counts)
        correl_nadir[i] <- sum(abs(pcm - scm))
      }
      
      # SHOULD WE SAVE AND RETURN THE SIMULATED VALUES?
      # ASR: We compute the mean simulated value and return it as an attribute
      #      because we want to explore the simulated values in the future.
      if (nadir[[2]]) {
        res <- list(strata = strata_nadir, correl = correl_nadir)
      } else {
        res <- list(strata = "strata_nadir", correl = "correl_nadir")
      }
      a <- attributes(res)
      a$strata_nadir <- mean(strata_nadir) / scale$max
      a$correl_nadir <- mean(correl_nadir) / scale$max
      attributes(res) <- a
      
      # ASR: Other options are not available yet.
    } else {
      if (!is.null(nadir[[3]])) {
        message("sorry but you cannot set the nadir point")
      } else {
        if (!is.null(nadir[[4]])) {
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
  function (sm, n.cov, strata, pcm, scm, nadir, weights, n.pts, utopia, scale) {
    
    counts <- lapply(1:n.cov, function (i)
      hist(sm[, i], strata[[1]][[i]], plot = FALSE)$counts)
    
    # ASR: We use the proportion of points per sampling strata
    #counts <- lapply(1:n.cov, function(i) counts[[i]] / n.pts * 100)
    counts <- lapply(1:n.cov, function(i) counts[[i]] / n.pts)
    counts <- sapply(1:n.cov, function (i) 
      sum(abs(counts[[i]] - strata[[2]][[i]])))
        
    # scale the objective function values
    if(scale$type == "upper-lower") {
      obj_cont <- sum(counts) - utopia$strata / 
        attr(nadir, "strata") - utopia$strata  
      obj_cor <- sum(abs(pcm - scm)) - utopia$correl / 
        attr(nadir, "correl") - utopia$correl
    } else if (scale$type == "upper") {
      obj_cont <- sum(counts) / attr(nadir, "strata")
      obj_cor <- sum(abs(pcm - scm)) / attr(nadir, "correl")
    }
    
    # aggregate the objective function values
    obj_cont <- obj_cont * weights[[1]]
    obj_cor <- obj_cor * weights[[2]]
    res <- obj_cont + obj_cor
    return (res)
  }
# INTERNAL FUNCTION - NADIR FOR FACTOR COVARIATES ##############################
.facNadir <-
  function (nadir, candi, n.candi, n.pts, n.cov, covars, pop.prop, pcm, scale) {
    
    if (!is.null(nadir[[1]]) && !is.null(nadir[[2]])) { # Simulate the nadir
      m <- paste("simulating ", nadir[[1]], " nadir values...", sep = "")
      message(m)
      
      # SET VARIABLES
      strata_nadir <- vector()
      correl_nadir <- vector()
      
      # BEGIN THE SIMULATION
      for (i in 1:nadir[[1]]) {
        pts <- sample(1:n.candi, n.pts)
        sm <- covars[pts, ]
        scm <- pedometrics::cramer(sm)
        
        # ASR: We multiply the proportion by 100 for numerical stability
        samp_prop <- lapply(sm, function(x) table(x) / n.pts)
        #samp_prop <- lapply(sm, function(x) table(x) / n.pts * 100)
        samp_prop <- sapply(1:n.cov, function (i)
          sum(abs(samp_prop[[i]] - pop.prop[[i]])))
        strata_nadir[i] <- sum(samp_prop)
        correl_nadir[i] <- sum(abs(pcm - scm))
      }
      
      # SHOULD WE SAVE AND RETURN THE SIMULATED VALUES?
      # ASR: We compute the mean simulated value and return it as an attribute
      #      because we want to explore the simulated values in the future.
      if (nadir[[2]]) {
        res <- list(strata = strata_nadir, correl = correl_nadir)
      } else {
        res <- list(strata = "strata_nadir", correl = "correl_nadir")
      }
      a <- attributes(res)
      a$strata_nadir <- mean(strata_nadir) / scale$max
      a$correl_nadir <- mean(correl_nadir) / scale$max
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
# INTERNAL FUNCTION - CRITERION FOR FACTOR COVARIATES ##########################
.objFac <-
  function (sm, pop.prop, nadir, weights, pcm, scm, n.pts, n.cov, utopia, 
            scale) {
    
    # ASR: We multiply the proportions by 100 for numerical stability
    samp_prop <- lapply(sm, function(x) table(x) / n.pts)
    #samp_prop <- lapply(sm, function(x) table(x) / n.pts * 100)
    samp_prop <- sapply(1:n.cov, function (i)
      sum(abs(samp_prop[[i]] - pop.prop[[i]])))
    
    # scale the objective function values
    if(scale$type == "upper-lower") {
      obj_cat <- sum(samp_prop) - utopia$strata / 
        attr(nadir, "strata") - utopia$strata
      obj_cor <- sum(abs(pcm - scm)) - utopia$correl / 
        attr(nadir, "correl") - utopia$correl
    } else if (scale$type == "upper") {
      obj_cat <- sum(samp_prop) / attr(nadir, "strata")
      obj_cor <- sum(abs(pcm - scm)) / attr(nadir, "correl")
    }
    
    # aggregate the objective function values
    obj_cat <- obj_cat * weights[[1]]
    obj_cor <- obj_cor * weights[[2]]
    res <- obj_cat + obj_cor
    return (res)
  }
# INTERNAL FUNCTION - UTOPIA POINT FOR FACTOR COVARIATES #######################
.facUtopia <-
  function (utopia) {
    if (!is.null(unlist(utopia$user))) {
      utopia <- list(correl = utopia$user$correl, strata = utopia$user$strata)
    } else {
      message("sorry but the utopia point cannot be calculated")
    }
  }
# INTERNAL FUNCTION - UTOPIA POINT FOR NUMERIC COVARIATES ######################
.numUtopia <-
  function (utopia) {
    if (!is.null(unlist(utopia$user))) {
      utopia <- list(correl = utopia$user$correl, strata = utopia$user$strata)
    } else {
      message("sorry but the utopia point cannot be calculated")
    }
  }
# CALCULATE OBJECTIVE FUNCTION VALUE ###########################################
#' @rdname optimACDC
#' @export
objACDC <-
  function (points, candi, covars, strata.type = "area", 
            weights = list(correl = 0.5, strata = 0.5), use.coords = FALSE, 
            utopia = list(user = list(correl = NULL, strata = NULL),
                          abs = NULL),
            nadir = list(sim = 1000, save.sim = TRUE, user = NULL, abs = NULL),
            scale = list(type = "upper-lower", max = 100)) {
    
    if (!is.data.frame(covars)) covars <- as.data.frame(covars)
    # Check arguments
    check <- .optimACDCcheck(candi = candi, covars = covars, nadir = nadir,
                             covars.type = covars.type, weights = weights, 
                             use.coords = use.coords, strata.type = strata.type,
                             utopia = utopia, scale = scale)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare sample points
    n_candi <- nrow(candi)
    points <- .spsannPoints(points = points, candi = candi, n.candi = n_candi)
    n_pts <- nrow(points)
    
    # Prepare covariates (covars) and create the starting sample matrix (sm)
    covars.type <- ifelse(pedometrics::is.any.factor(covars), "factor",
                          "numeric")
    covars <- .covarsACDC(covars = covars, covars.type = covars.type, 
                          use.coords = use.coords, candi = candi, n.pts = n_pts,
                          strata.type = strata.type)
    n_cov <- ncol(covars)
    sm <- covars[points[, 1], ]
    
    # Calculate the energy state
    if (covars.type == "numeric") { # Numeric covariates
      pcm <- cor(covars, use = "complete.obs")
      strata <- .numStrata(n.pts = n_pts, covars = covars, 
                           strata.type = strata.type)
      nadir <- .numNadir(n.pts = n_pts, n.cov = n_cov, n.candi = n_candi, 
                         pcm = pcm, nadir = nadir, candi = candi, 
                         covars = covars, strata = strata, scale = scale)
      utopia <- .numUtopia(utopia = utopia)
      scm <- cor(sm, use = "complete.obs")
      energy <- .objNum(sm = sm, n.cov = n_cov, strata = strata, pcm = pcm, 
                         scm = scm, nadir = nadir, weights = weights, 
                         n.pts = n_pts, utopia = utopia, scale = scale)
      
    } else { # Factor covariates
      if (covars.type == "factor") {
        pcm <- pedometrics::cramer(covars)
        #pop_prop <- lapply(covars, function(x) table(x) / nrow(covars) * 100)
        pop_prop <- lapply(covars, function(x) table(x) / nrow(covars))
        nadir <- .facNadir(nadir = nadir, candi = candi, n.candi = n_candi,
                           n.pts = n_pts, n.cov = n_cov, covars = covars, 
                           pop.prop = pop_prop, pcm = pcm, scale = scale)
        utopia <- .facUtopia(utopia = utopia)
        scm <- pedometrics::cramer(sm)
        energy <- .objFac(sm = sm, pop.prop = pop_prop, nadir = nadir, 
                           weights = weights, pcm = pcm, scm = scm,
                           n.pts = n_pts, n.cov = n_cov, utopia = utopia,
                           scale = scale)
      }
    }
    return (energy)
  }
# INTERNAL FUNCTION - USE THE COORDINATES AS COVARIATES ########################
.useCoords <-
  function (covars.type, candi, n.pts, strata.type) {
    if (covars.type == "factor") {
      coords <- data.frame(candi[, 2:3])
      breaks <- .numStrata(n.pts = n.pts, covars = coords, 
                           strata.type = strata.type)[[1]]
      coords <- cont2cat(x = coords, breaks = breaks)
      covars <- data.frame(covars, coords)
    } else {
      covars <- data.frame(covars, candi[, 2:3])
    }
    return (covars)
  }
# INTERNAL FUNCTION - PREPARE THE COVARIATES ###################################
.covarsACDC <-
  function (covars, covars.type, use.coords, candi, n.pts, strata.type) {
    
    # Factor covariates
    if (covars.type == "factor") {
      if (use.coords) {
        covars <- data.frame(covars, candi[, 2:3])
      }
      # Convert numeric covariates to factor covariates
      if (!pedometrics::is.all.factor(covars)) {
        i <- which(sapply(covars, is.factor) == FALSE)
        mes <- paste("converting ", length(i), 
                     " numeric covariates to factor covariates", sep = "")
        message(mes)
        num_covars <- data.frame(covars[, i])
        breaks <- .numStrata(n.pts = n.pts, covars = num_covars, 
                             strata.type = strata.type)[[1]]
        num_covars <- cont2cat(x = num_covars, breaks = breaks)
        covars[, i] <- num_covars
      }
      
      # Numeric covariates
    } else {
      if (use.coords) {
        covars <- data.frame(covars, candi[, 2:3])
      }
    }
    return (covars)
  }
