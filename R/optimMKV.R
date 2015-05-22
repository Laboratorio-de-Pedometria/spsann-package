#' Optimization of sample configurations for spatial interpolation
#'
#' Optimize a sample configuration for spatial interpolation with a known linear 
#' model. A criterion is defined so that the sample configuration minimizes the
#' mean/maximum kriging variance.
#'
#' @template spJitter_doc
#' @template spSANN_doc
#' 
#' @param covars Data frame or matrix with the covariates in the columns.
#' 
#' @param equation Formula string that defines the dependent variable \code{z}
#' as a linear model of independent variables. Defaults to 
#' \code{equation = z ~ 1}. See the argument \code{formula} in the function
#' \code{\link[gstat]{krige}} for more information.
#'
#' @param model Object of class "variogramModel". See the argument 
#' \code{model} in the function \code{\link[gstat]{krige}} for more information.
#'
#' @param krige.stat Character value defining the statistic that should be used
#' to summarize the kriging variance. Available options are \code{"mean"} and
#' \code{"max"} for the mean and maximum kriging variance, respectively.
#' Defaults to \code{krige.stat = "mean"}.
#'
#' @return
#' \code{optimMKV} returns a matrix: the optimized sample configuration with
#' the evolution of the energy state during the optimization as an attribute.
#'
#' \code{objMKV} returns a numeric value depending on the choice of 
#' \code{krige.stat}. If \code{krige.stat = "mean"}, the mean kriging variance.
#' If \code{krige.stat = "max"}, the maximum kriging variance.
#'
#' @references
#' Brus, D. J. & Heuvelink, G. B. M. Optimization of sample patterns for
#' universal kriging of environmental variables. \emph{Geoderma}. v. 138, 
#' p. 86-95, 2007.
#' 
#' Heuvelink, G. B. M.; Brus, D. J. & de Gruijter, J. J. Optimization of sample
#' configurations for digital mapping of soil properties with universal kriging.
#' In: Lagacherie, P.; McBratney, A. & Voltz, M. (Eds.) \emph{Digital soil
#' mapping - an introductory perspective}. Elsevier, v. 31, p. 137-151, 2006.
#' 
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @aliases optimMKV
#' @keywords spatial optimize
#' @concept simulated annealing
#' @importFrom plyr is.formula
#' @export
#' 
# FUNCTION - MAIN ##############################################################
optimMKV <-
  function (points, candi, covars, equation = z ~ 1, model, krige.stat = "mean",
            x.max, x.min, y.max, y.min, iterations = 10000,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            boundary, progress = TRUE, verbose = TRUE, greedy = FALSE) {
    
    if (!missing(covars)) {
      if (!is.data.frame(covars)) covars <- as.data.frame(covars) 
    }    
    
    # Check arguments
    check <- .spSANNcheck(points = points, candi = candi, x.max = x.max, 
                          x.min = x.min, y.max = y.max, y.min = y.min, 
                          iterations = iterations, acceptance = acceptance,
                          stopping = stopping, plotit = plotit, 
                          boundary = boundary, progress = progress, 
                          verbose = verbose)
    if (!is.null(check)) stop (check, call. = FALSE)
    check <- .optimMKVcheck(covars = covars, equation = equation, model = model,
                            krige.stat = krige.stat, candi = candi)
    if (!is.null(check)) stop (check, call. = FALSE)
    
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
    MOOP <- FALSE
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
      #new_row <- covars[new_conf[wp, 1], ]
      new_row <- cbind(1, covars[new_conf[wp, 1], ])
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
                    x_max0 = x_max0, x.max = x.max, best.k = best_k, 
                    best.energy = best_energy, MOOP = FALSE)
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
    res <- .spSANNout(new_conf, energy0, energies, time0, MOOP = FALSE)
    return (res)
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
.optimMKVcheck <-
  function (covars, equation, model, krige.stat, candi) {
    
    # covars
    if (!missing(covars)) {
      if (nrow(candi) != nrow(covars)) {
        res <-
          paste("'candi' and 'covars' must have the same number of rows")
        return (res)
      }
    }
    
    # equation
    bb <- !is.formula(equation)
    cc <- all.vars(equation)[1] != "z"
    if (bb || cc) {
      res <-
        paste("'equation' must be a formula with dependent variable 'z'")
      return (res)
    }
    
    # model
    aa <- missing(model)
    bb <- class(model)[1] != "variogramModel"
    if (aa || bb) {
      res <- paste("'model' must be of class 'variogramModel'")
      return (res)
    }
    
    # krige.stat
    
    aa <- match(krige.stat, c("mean", "max"))
    if (is.na(aa)) {
      res <- paste("'krige.stat = ", krige.stat, "' is not supported",
                   sep = "")
      return (res)
    }
  }
# FUNCTION - CANCLULATE THE OBJECTIVE FUNCTION VALUE ###########################
#' @export
#' @rdname optimMKV
objMKV <-
  function (points, candi, covars, equation, model, krige.stat = "mean") {
    
    if (!missing(covars)) {
      if (!is.data.frame(covars)) covars <- as.data.frame(covars) 
    }    
    
    # Check arguments
    #check <- .spSANNcheck(points, candi)
    #if (!is.null(check)) stop (check, call. = FALSE)
    
    #check <- .optimMKVcheck()
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
