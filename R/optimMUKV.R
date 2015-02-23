#' Optimization of sample patterns with a known model
#'
#' Optimize a sample pattern with a known model.
#'
#' @template spJitter_doc
#' @template spSANN_doc
#' 
#' @param lags Integer. The number of lag distance classes. Alternatively, a
#' vector of numeric values with the lower and upper limits of each lag
#' distance class. The lowest value must be larger than zero, e.g. 0.0001.
#' Defaults to \code{lags = 7}.
#' 
#' @param lags.type Character. The type of lag distance classes. Available
#' options are \code{"equidistant"} and \code{"exponential"}. Defaults to
#' \code{lags.type = "exponential"}. See \sQuote{Details} for more information.
#'
#' @param lags.base Numeric. Base of the exponential expression used to
#' create the exponential lag distance classes. Defaults to
#' \code{lags.base = 2}. See \sQuote{Details} for more information.
#'
#' @param cutoff Numeric value. The maximum distance up to which lag distance
#' classes are created. Used only when \code{lags} is an integer. 
#' See \sQuote{Details} for more information.
#'
#' @param criterion Character value. The measure to be used to describe the
#' energy state of the current system configuration. Available options are
#' \code{"minimum"} and \code{"distribution"}. Defaults to
#' \code{objective = "minimum"}. See \sQuote{Details} for more information.
#'
#' @param pre.distri Numeric vector. The pre-specified distribution of points
#' or point-pair with which the observed counts of points or point-pairs per
#' lag distance class is compared. Used only when
#' \code{criterion = "distribution"}. Defaults to a uniform distribution. See
#' \sQuote{Details} for more information.
#'
#' @return
#' \code{optimPPL} returns a matrix: the optimized sample pattern with
#' the evolution of the energy state during the optimization as an attribute.
#'
#' \code{pointsPerLag} and \code{pairsPerLag} return a data.frame with three
#' columns: a) the lower and b) upper limits of each lag, and c) the number of
#' points or point-pairs per lag.
#'
#' \code{objPoints} and \code{objPairs} return a numeric value depending on the
#' choice of \code{criterion}. If \code{criterion = "distribution"}, the sum of
#' the differences between the pre-specified and observed distribution of counts
#' of points or point-pairs per lag. If \code{criterion = "minimum"}, the
#' inverse of the minimum count of points or point pairs over all lags
#' multiplied by a constant.
#'
#' @references
#' Bresler, E.; Green, R. E. \emph{Soil parameters and sampling scheme for
#' characterizing soil hydraulic properties of a watershed}. Honolulu:
#' University of Hawaii at Manoa, p. 42, 1982.
#'
#' Marler, R. T.; Arora, J. S. Function-transformation methods for
#' multi-objective optimization. \emph{Engineering Optimization}. v. 37, p.
#' 551-570, 2005.
#'
#' Russo, D. Design of an optimal sampling network for estimating the variogram.
#' \emph{Soil Science Society of America Journal}. v. 48, p. 708-716, 1984.
#'
#' Truong, P. N.; Heuvelink, G. B. M.; Gosling, J. P. Web-based tool for expert
#' elicitation of the variogram. \emph{Computers and Geosciences}. v. 51, p.
#' 390-399, 2013.
#'
#' Warrick, A. W.; Myers, D. E. Optimization of sampling locations for variogram
#' calculations. \emph{Water Resources Research}. v. 23, p. 496-500, 1987.
#'
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @aliases optimMUKV
#' @keywords spatial optimize
#' @concept simulated annealing
#' @export
#' @examples
#' require(sp)
#' data(meuse)
#' meuse <- as.matrix(meuse[, 1:2])
#' meuse <- matrix(cbind(c(1:dim(meuse)[1]), meuse), ncol = 3)
#' pointsPerLag(meuse, cutoff = 1000)
#' objPoints(meuse, cutoff = 1000)
# FUNCTION - MAIN ##############################################################
optimMUKV <-
  function (points, candi, covars, equation = z ~ 1, model, krige.stat = "mean",
            x.max, x.min, y.max, y.min, iterations = 10000,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            boundary, progress = TRUE, verbose = TRUE) {
    
    if (!missing(covars)) {
      if (!is.data.frame(covars)) covars <- as.data.frame(covars) 
    }    
    
    # Check arguments
    # http://www.r-bloggers.com/a-warning-about-warning/
    check <- .spSANNcheck(points, candi, x.max, x.min, y.max, y.min,
                          iterations, acceptance, stopping, plotit, boundary,
                          progress, verbose)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    #check <- .optimMUKVcheck()
    #if (!is.null(check)) stop (check, call. = FALSE)
    
    if (plotit) {
      par0 <- par()
      on.exit(suppressWarnings(par(par0)))
    }
    
    # Prepare points
    n_candi <- nrow(candi)
    points <- .spsannPoints(points = points, candi = candi, n.candi = n_candi)
    n_pts <- nrow(points)
    conf0 <- points
    old_conf <- conf0
    
    # Prepare prediction grid (pg) with covars
    if (terms(equation)[[3]] == 1) {
      covars <- data.frame(candi[, 2:3])
      colnames(covars) <- c("x", "y")
    } else {
      covars <- data.frame(candi[, 2:3], covars[, all.vars(equation)[-1]])
      colnames(covars) <- c("x", "y", all.vars(equation)[-1])
    }
    
    # Prepare starting sample matrix (sm)
    z <- rep(1, n_pts)
    if (terms(equation)[[3]] == 1) {
      sm <- data.frame(z, points[, 2:3])
      colnames(sm) <- c("z", "x", "y")
    } else {
      sm <- data.frame(z, points[, 2:3],
                       covars[points[, 1], all.vars(equation)[-1]])
      colnames(sm) <- c("z", "x", "y", all.vars(equation)[-1])
    }
    
    # Initial energy state
    energy0 <- gstat::krige(formula = equation, locations = ~ x + y, data = sm,
                            newdata = covars, model = model)$var1.var
    if (krige.stat == "mean") {
      energy0 <- mean(energy0)
    } else {
      if (krige.stat == "max") {
        energy0 <- max(energy0)
      }
    }
    
    # other settings for the simulated annealing algorithm
    old_sm <- sm
    new_sm <- sm
    best_sm <- sm
    count <- 0
    old_energy <- energy0
    best_energy <- Inf
    energies <- vector()
    accept_probs <- vector()
    x_max0 <- x.max
    y_max0 <- y.max
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()
    
    # begin the main loop
    for (k in 1:iterations) {
      
      # jitter one of the points and update x.max and y.max
      wp <- sample(1:n_pts, 1)
      new_conf <- spJitterFinite(points = old_conf, candi = candi,
                                 x.max = x.max, x.min = x.min, y.max = y.max,
                                 y.min = y.min, which.point = wp)
      x.max <- x_max0 - (k / iterations) * (x_max0 - x.min)
      y.max <- y_max0 - (k / iterations) * (y_max0 - y.min)
      
      # Update sample matrix and energy state
      new_row <- covars[new_conf[wp, 1], ]
      new_sm[wp, ] <- new_row
      new_energy <- gstat::krige(formula = equation, locations = ~ x + y, 
                                 data = new_sm, newdata = covars, 
                                 model = model, debug.level = 0)$var1.var
      if (krige.stat == "mean") {
        new_energy <- mean(new_energy)
      } else {
        if (krige.stat == "max") {
          new_energy <- max(new_energy)
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
      if (plotit && any(round(seq(1, iterations, 10)) == k)) {
        .spSANNplot(energy0 = energy0, energies = energies, k = k, 
                    acceptance = acceptance, accept_probs = accept_probs,
                    boundary = boundary, new_conf = new_conf[, 2:3], 
                    conf0 = conf0[, 2:3], y_max0 = y_max0, y.max = y.max, 
                    x_max0 = x_max0, x.max = x.max)
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
.optimmukvcheck <-
  function () {
    
  }
# FUNCTION - CANCLULATE THE OBJECTIVE FUNCTION VALUE ###########################
objMUKV <-
  function (points, candi, covars, equation, model, krige.stat = "mean") {
    
    if (!missing(covars)) {
      if (!is.data.frame(covars)) covars <- as.data.frame(covars) 
    }    
    
    # Check arguments
    #check <- .spSANNcheck(points, candi)
    #if (!is.null(check)) stop (check, call. = FALSE)
    
    #check <- .optimMUKVcheck()
    #if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points
    n_candi <- nrow(candi)
    points <- .spsannPoints(points = points, candi = candi, n.candi = n_candi)
    n_pts <- nrow(points)
    
    # Prepare prediction grid (pg) with covars
    if (terms(equation)[[3]] == 1) {
      covars <- data.frame(candi[, 2:3])
      colnames(covars) <- c("x", "y")
    } else {
      covars <- data.frame(candi[, 2:3], covars[, all.vars(equation)[-1]])
      colnames(covars) <- c("x", "y", all.vars(equation)[-1])
    }
    
    # Prepare starting sample matrix (sm)
    z <- rep(1, n_pts)
    if (terms(equation)[[3]] == 1) {
      sm <- data.frame(z, points[, 2:3])
      colnames(sm) <- c("z", "x", "y")
    } else {
      sm <- data.frame(z, points[, 2:3],
                       covars[points[, 1], all.vars(equation)[-1]])
      colnames(sm) <- c("z", "x", "y", all.vars(equation)[-1])
    }
    
    # Initial energy state
    res <- gstat::krige(formula = equation, locations = ~ x + y, data = sm,
                        newdata = covars, model = model)$var1.var
    if (krige.stat == "mean") {
      res <- mean(res)
    } else {
      if (krige.stat == "max") {
        res <- max(res)
      }
    }
    return (res)
  }
