#' Optimization of sample configurations for spatial interpolation
#'
#' Optimize a sample configuration for spatial interpolation with a known linear
#' model. A criterion is defined so that the sample configuration minimizes the
#' mean/maximum kriging variance (\bold{MKV}).
#'
#' @template spJitter_doc
#' @template spSANN_doc
#' @template MOOP_doc
#' 
#' @param covars Data frame or matrix with the covariates in the columns.
#' 
#' @param eqn Formula string that defines the dependent variable \code{z}
#' as a linear model of the independent variables contained in \code{covars}. 
#' Defaults to \code{eqn = z ~ 1}, that is, ordinary kriging. See the 
#' argument \code{formula} in the function \code{\link[gstat]{krige}} for 
#' more information.
#'
#' @param vgm Object of class "variogramModel". See the argument 
#' \code{model} in the function \code{\link[gstat]{krige}} for more 
#' information.
#'
#' @param krige.stat Character value defining the statistic that should be used
#' to summarize the kriging variance. Available options are \code{"mean"} and
#' \code{"max"} for the mean and maximum kriging variance, respectively.
#' Defaults to \code{krige.stat = "mean"}.
#' 
#' @param ... further arguments passed to \code{\link[gstat]{krige}}.
#'
#' @return
#' \code{optimMKV} returns a matrix: the optimized sample configuration.
#'
#' \code{objMKV} returns a numeric value: the energy state of the sample
#' configuration - the objective function value.
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
#' @aliases optimMKV objMKV
#' @concept spatial interpolation
#' @export
#' @examples
#' require(sp)
#' require(gstat)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' covars <- as.data.frame(meuse.grid)
#' vgm <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
#' set.seed(2001)
#' res <- optimMKV(points = 100, candi = candi, covars = covars, maxdist = 500,
#'                 eqn = z ~ dist, vgm = vgm)
#' tail(attr(res, "energy"), 1) # 11.9878
#' objMKV(points = res, candi = candi, covars = covars, eqn = z ~ dist, 
#'        vgm = vgm, maxdist = 500)
# FUNCTION - MAIN ##############################################################
optimMKV <-
  function (
    # MKV
    covars, eqn = z ~ 1, vgm, krige.stat = "mean", ...,
    # SPSANN
    points, candi, iterations = 100, x.max, x.min, y.max, y.min,
    acceptance = list(initial = 0.99, cooling = iterations / 10),
    stopping = list(max.count = iterations / 10), plotit = FALSE, track = FALSE,
    boundary, progress = TRUE, verbose = FALSE, greedy = FALSE,
    # MOOP
    weights = NULL, nadir = NULL, utopia = NULL) {
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Check other arguments
    check <- .checkMKV(covars = covars, eqn = eqn, vgm = vgm, 
                       krige.stat = krige.stat, candi = candi)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Prepare prediction grid (covars) and starting sample matrix (sm)
    covars <- .covarsMKV(eqn = eqn, candi = candi, covars = covars)
    sm <- .smMKV(n_pts = n_pts, eqn = eqn, pts = points, covars = covars)
    
    # Initial energy state
    energy0 <- .objMKV(eqn = eqn, sm = sm, covars = covars, vgm = vgm, 
                       krige.stat = krige.stat, ...)
    
    # other settings for the simulated annealing algorithm
    old_sm <- sm
    new_sm <- sm
    best_sm <- sm
    count <- 0
    old_energy <- energy0
    best_energy <- Inf
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()
    
    # begin the main loop
    for (k in 1:iterations) {
      
      # Plotting and jittering
      eval(.plot_and_jitter())
      
      # Update sample matrix and energy state
      new_sm[wp, ] <- cbind(1, covars[new_conf[wp, 1], ])
      new_energy <- .objMKV(eqn = eqn, sm = new_sm, covars = covars, vgm = vgm, 
                            krige.stat = krige.stat, debug.level = 0, ...)
      
      # Evaluate the new system configuration
      if (greedy) {
        random_prob <- 1
      } else {
        random_prob <- runif(1)
      }
      actual_prob <- acceptance[[1]] * exp(-k / acceptance[[2]])
      if (track) accept_probs[k] <- actual_prob
      if (new_energy <= old_energy) {
        old_conf <- new_conf
        old_energy <- new_energy
        count <- 0
        old_sm <- new_sm
      } else {
        if (new_energy > old_energy & random_prob <= actual_prob) {
          old_conf <- new_conf
          old_energy <- new_energy
          count <- count + 1
          old_sm <- new_sm
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... p = ",
                random_prob, "\n")
          }
        } else {
          new_energy <- old_energy
          new_conf <- old_conf
          count <- count + 1
          new_sm <- old_sm
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... stops at",
                stopping[[1]], "\n")
          }
        }
      }
      
      # Best energy state
      if (track) energies[k] <- new_energy
      if (new_energy < best_energy / 1.0000001) {
        best_k <- k
        best_conf <- new_conf
        best_energy <- new_energy
        best_old_energy <- old_energy
        old_conf <- old_conf
        best_sm <- new_sm
        best_old_sm <- old_sm
      }
      
      # Freezing parameters
      if (count == stopping[[1]]) {
        if (new_energy > best_energy * 1.000001) {
          old_conf <- old_conf
          new_conf <- best_conf
          old_energy <- best_old_energy
          new_energy <- best_energy
          count <- 0
          new_sm <- best_sm
          old_sm <- best_old_sm
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
    
    # Prepare output
    eval(.prepare_output())
  }
# INTERNAL FUNCTION - CALCULATE THE ENERGY STATE VALUE #########################
# eqn: equation
# sm: sample matrix
# covars: covariates
# vgm: variogram model
# krige.stat: statistic to summarize the kriging variance
.objMKV <-
  function (eqn, sm, covars, vgm, krige.stat, ...) {
    
    res <- gstat::krige(formula = eqn, locations = ~ x + y, data = sm,
                        newdata = covars, model = vgm, ...)$var1.var
    
    # Calculate the energy state value
    if (krige.stat == "mean") { # Mean kriging variance
      res <- mean(res)
      
    } else { # Maximum kriging variance
      res <- max(res)
    }
    
    # Output
    return (res)
  }
# INTERNAL FUNCTION - PREPARE STARTING SAMPLE MATRIX ###########################
# n_pts: number of points
# eqn: equation
# pts: points
# covars: covariates
.smMKV <-
  function (n_pts, eqn, pts, covars) {
    
    z <- rep(1, n_pts)
    if (terms(eqn)[[3]] == 1) { # Simple and ordinary kriging
      sm <- data.frame(z, pts[, 2:3])
      colnames(sm) <- c("z", "x", "y")
      
    } else { # Universal kriging
      sm <- data.frame(z, pts[, 2:3],
                       covars[pts[, 1], all.vars(eqn)[-1]])
      colnames(sm) <- c("z", "x", "y", all.vars(eqn)[-1])
    }
    
    # Output
    return (sm)
  }
# INTERNAL FUNCTION - PREPARE THE PREDICITON GRID ##############################
# eqn: equation
# candi: candidate locations
# covars: covariates
.covarsMKV <-
  function (eqn, candi, covars) {
    
    if (terms(eqn)[[3]] == 1) { # Simple and ordinary kriging
      covars <- data.frame(candi[, 2:3])
      colnames(covars) <- c("x", "y")
      
    } else { # Universal kriging
      covars <- as.data.frame(covars)
      covars <- data.frame(candi[, 2:3], covars[, all.vars(eqn)[-1]])
      colnames(covars) <- c("x", "y", all.vars(eqn)[-1])
    }
    
    # Output
    return (covars)
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
.checkMKV <-
  function (covars, eqn, vgm, krige.stat, candi) {
    
    # covars
    if (!missing(covars)) {
      if (is.vector(covars)) {
        if (nrow(candi) != length(covars)) {
          res <- "'candi' and 'covars' must have the same number of rows"
          return (res)
        }
      } else {
        if (nrow(candi) != nrow(covars)) {
          res <- "'candi' and 'covars' must have the same number of rows"
          return (res)
        }
      }
    }
    
    # eqn
    bb <- !inherits(eqn, "formula")
    cc <- all.vars(eqn)[1] != "z"
    if (bb || cc) {
      res <- "'eqn' must be a formula with dependent variable 'z'"
      return (res)
    }
    
    # vgm
    aa <- missing(vgm)
    bb <- class(vgm)[1] != "variogramModel"
    if (aa || bb) {
      res <- "'vgm' must be of class 'variogramModel'"
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
  function (points, candi, covars, eqn, vgm, krige.stat = "mean", ...) {
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare prediction grid with covars
    covars <- .covarsMKV(eqn = eqn, candi = candi, covars = covars)
    sm <- .smMKV(n_pts = n_pts, eqn = eqn, pts = points, covars = covars)
    
    # Energy state
    res <- .objMKV(eqn = eqn, sm = sm, covars = covars, vgm = vgm, 
                   krige.stat = krige.stat, ...)
    
    # Output
    return (res)
  }
