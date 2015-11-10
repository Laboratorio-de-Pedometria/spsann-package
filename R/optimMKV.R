#' Optimization of sample configurations for spatial interpolation
#'
#' Optimize a sample configuration for spatial interpolation with a known linear
#' model. A criterion is defined so that the sample configuration minimizes the
#' mean/maximum kriging variance (\bold{MKV}).
#'
#' @inheritParams spJitterFinite
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
#' @note
#' This function is based on the method originally proposed by Heuvelink, Brus 
#' and de Gruijter (2006) and implemented in the R-package 
#' \pkg{intamapInteractive} by Edzer Pebesma and Jon Skoien.
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
#' schedule <- scheduleSPSANN(initial.temperature = 0.5, chains = 1)
#' \dontrun{
#' # This example takes more than 5 seconds to run!
#' set.seed(1984)
#' res <- optimMKV(points = 100, candi = candi, covars = covars, maxdist = 600,
#'                 eqn = z ~ dist, vgm = vgm, schedule = schedule)
#' objSPSANN(res) -
#'   objMKV(points = res, candi = candi, covars = covars, eqn = z ~ dist, 
#'          vgm = vgm, maxdist = 600)
#' }
# FUNCTION - MAIN ##############################################################
optimMKV <-
  function (points, candi,
            # MKV
            covars, eqn = z ~ 1, vgm, krige.stat = "mean", ...,
            # SPSANN
            schedule = scheduleSPSANN(), plotit = FALSE, track = FALSE,
            boundary, progress = TRUE, verbose = FALSE,
            # MOOP
            weights = NULL, nadir = NULL, utopia = NULL) {
    
    # Check suggests
    pkg <- c("gstat")
    eval(.check_suggests())
    
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
    
    # Other settings for the simulated annealing algorithm
    old_sm <- sm
    new_sm <- sm
    best_sm <- sm
    old_energy <- energy0
    best_energy <- Inf
    actual_temp <- schedule$initial.temperature
    k <- 0 # count the number of jitters
    if (progress) {
      max <- n_pts * schedule$chains * schedule$chain.length
      pb <- utils::txtProgressBar(min = 1, max = max, style = 3)
    }
    time0 <- proc.time()
    
    # Initiate the annealing schedule
    for (i in 1:schedule$chains) {
      n_accept <- 0
      
      for (j in 1:schedule$chain.length) { # Initiate one chain
        
        for (wp in 1:n_pts) { # Initiate loop through points
          k <- k + 1
          
          # Plotting and jittering
          eval(.plot_and_jitter())
          
          # Update sample matrix and energy state
          new_sm[wp, ] <- cbind(1, covars[new_conf[wp, 1], ])
          new_energy <- .objMKV(eqn = eqn, sm = new_sm, covars = covars, 
                                vgm = vgm, krige.stat = krige.stat, 
                                debug.level = 0, ...)
          
          # Evaluate the new system configuration
          accept <- min(1, exp((old_energy - new_energy) / actual_temp))
          accept <- floor(rbinom(n = 1, size = 1, prob = accept))
          if (accept) {
            old_conf <- new_conf
            old_energy <- new_energy
            old_sm <- new_sm
            n_accept <- n_accept + 1
          } else {
            new_energy <- old_energy
            new_conf <- old_conf
            new_sm <- old_sm
          }
          if (track) energies[k] <- new_energy
          
          # Record best energy state
          if (new_energy < best_energy / 1.0000001) {
            best_k <- k
            best_conf <- new_conf
            best_energy <- new_energy
            best_old_energy <- old_energy
            old_conf <- old_conf
            best_sm <- new_sm
            best_old_sm <- old_sm
          }
          
          
          if (progress) utils::setTxtProgressBar(pb, k)
        } # End loop through points
        
      } # End the chain
      
      # Check the proportion of accepted swaps in the first chain
      if (i == 1) {
        x <- round(n_accept / c(n_pts * schedule$chain.length), 2)
        if (x < schedule$initial.acceptance) {
          cat("\nlow temperature: only ", x," of acceptance in the 1st chain\n", 
              sep = "")
          break
        }
      }
      
      # Count the number of chains without any change in the objective function.
      # Restart with the previously best configuration if it exists.
      if (n_accept == 0) {
        no_change <- no_change + 1
        if (no_change > schedule$stopping) {
          if (new_energy > best_energy * 1.000001) {
            old_conf <- old_conf
            new_conf <- best_conf
            old_energy <- best_old_energy
            new_energy <- best_energy
            new_sm <- best_sm
            old_sm <- best_old_sm
            no_change <- 0
            cat("\nrestarting with previously best configuration\n")
          } else { 
            break 
          }
        }
        if (verbose) {
          cat("\n", no_change, "chain(s) with no improvement... stops at",
              schedule$stopping, "\n")
        }
      } else {
        no_change <-  0
      }
      
      # Update control parameters
      actual_temp <- actual_temp * schedule$temperature.decrease
      x.max <- x_max0 - (i / schedule$chains) * (x_max0 - x.min)
      y.max <- y_max0 - (i / schedule$chains) * (y_max0 - y.min)
      
    } # End the annealing schedule
    
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
    if (stats::terms(eqn)[[3]] == 1) { # Simple and ordinary kriging
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
    
    if (stats::terms(eqn)[[3]] == 1) { # Simple and ordinary kriging
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
  function (points, candi,
            # MKV
            covars, eqn = z ~ 1, vgm, krige.stat = "mean", ...) {
    
    # Check suggests
    pkg <- c("gstat")
    eval(.check_suggests())
    
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
