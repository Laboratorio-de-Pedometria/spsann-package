#' Optimization of sample configurations for spatial interpolation (II)
#'
#' Optimize a sample configuration for spatial interpolation with a 'known' linear mixed model, e.g. universal 
#' (external drift) kriging and regression-kriging with a linear regression model. A criterion is defined so
#' that the sample configuration minimizes the mean or maximum kriging prediciton error variance (\bold{MKV}).
#'
#' @template spSANN_doc
#' @template spJitter_doc
#' @template schedule_doc
#'
# @inheritParams optimMSSD
#'
#' @param covars Data frame or matrix with the covariates in the columns. The number of rows of `covars` must 
#' match exactly that of `candi` -- or `eval.grid`, in case a coarser evaluation grid is used.
#' 
#' @param eqn Formula string that defines the dependent variable `z` as a linear function of the independent 
#' variables (covariates) contained in `covars`. See the argument `formula` in the function 
#' `\link[gstat]{krige}` for more information.
#'
#' @param vgm Object of class `variogramModel`. See the argument `model` in the function `\link[gstat]{krige}`
#' for more information.
#'
#' @param krige.stat Character value defining the statistic that should be used to summarize the kriging 
#' prediction error variance. Available options are `"mean"` and `"max"` for the mean and maximum kriging  
#' prediction error variance, respectively. Defaults to `krige.stat = "mean"`.
#' 
#' @param ... further arguments passed to `\link[gstat]{krige}`. (Advanced users only!)
#'
#' @return
#' \code{optimMKV} returns an object of class \code{OptimizedSampleConfiguration}: the optimized sample
#' configuration with details about the optimization.
#'
#' \code{objMKV} returns a numeric value: the energy state of the sample configuration -- the objective
#' function value.
#'
#' @note
#' This function is based on the method originally proposed by Heuvelink, Brus and de Gruijter (2006) and 
#' implemented in the R-package \pkg{intamapInteractive} by Edzer Pebesma and Jon Skoien.
#'
#' @references
#' Brus, D. J.; Heuvelink, G. B. M. Optimization of sample patterns for universal kriging of environmental 
#' variables. \emph{Geoderma}. v. 138, p. 86-95, 2007.
#' 
#' Heuvelink, G. B. M.; Brus, D. J.; de Gruijter, J. J. Optimization of sample configurations for digital 
#' mapping of soil properties with universal kriging. In: Lagacherie, P.; McBratney, A. & Voltz, M. (Eds.)
#' \emph{Digital soil mapping - an introductory perspective}. Elsevier, v. 31, p. 137-151, 2006.
#' 
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @aliases optimMKV objMKV MKV
#' @concept spatial interpolation
#' @export
#' @examples
#' #####################################################################
#' # NOTE: The settings below are unlikely to meet your needs.         #
#' #####################################################################
#' \dontrun{
#' data(meuse.grid, package = "sp")
#' candi <- meuse.grid[1:1000, 1:2]
#' covars <- as.data.frame(meuse.grid)[1:1000, ]
#' vgm <- gstat::vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
#' schedule <- scheduleSPSANN(
#'   initial.temperature = 10, chains = 1, x.max = 1540, y.max = 2060, 
#'   x.min = 0,  y.min = 0, cellsize = 40)
#' set.seed(2001)
#' res <- optimMKV(
#'   points = 10, candi = candi, covars = covars, eqn = z ~ dist, 
#'   vgm = vgm, schedule = schedule)
#' data.frame(
#'   expected = 15.37137,
#'   objSPSANN = objSPSANN(res),
#'   objMKV = objMKV(
#'     points = res, candi = candi, covars = covars, eqn = z ~ dist, vgm = vgm)
#' )
#' }
# FUNCTION - MAIN #############################################################################################
optimMKV <-
  function (points, candi, 
            # eval.grid,
            # MKV
            covars, eqn, vgm, krige.stat = "mean", ...,
            # SPSANN
            schedule, plotit = FALSE, track = FALSE,
            boundary, progress = "txt", verbose = FALSE) {
    
    # Objective function name
    objective <- "MKV"
    
    # Check suggests
    pkg <- c("gstat")
    eval(.check_suggests())
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Check other arguments
    check <- .checkMKV(
      covars = covars, eqn = eqn, vgm = vgm, krige.stat = krige.stat, 
      # eval.grid = eval.grid,
      candi = candi)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Prepare prediction grid (covars) and starting sample matrix (sm)
    # if (!missing(eval.grid)) { # Use coarser prediction (evaluation) grid
      # covars <- .covarsMKV(eqn = eqn, pred.grid = eval.grid, covars = covars)
    # } else { # Use candi as prediction (evaluation) grid
      covars <- .covarsMKV(eqn = eqn, pred.grid = candi[, 2:3], covars = covars)
    # }
    sm <- .smMKV(n_pts = n_pts + n_fixed_pts, eqn = eqn, pts = points, covars = covars)
    
    # Initial energy state
    energy0 <- data.frame(
      obj = .objMKV(eqn = eqn, sm = sm, covars = covars, vgm = vgm, krige.stat = krige.stat, k = 0, ...))
    
    # Other settings for the simulated annealing algorithm
    old_sm <- sm
    new_sm <- sm
    best_sm <- sm
    old_energy <- energy0
    best_energy <- data.frame(obj = Inf)
    actual_temp <- schedule$initial.temperature
    k <- 0 # count the number of jitters
    
    # Set progress bar
    eval(.set_progress())

    # Initiate the annealing schedule
    for (i in 1:schedule$chains) {
      n_accept <- 0
      
      for (j in 1:schedule$chain.length) { # Initiate one chain
        
        for (wp in 1:n_pts) { # Initiate loop through points
          k <- k + 1

          # Plotting and jittering
          eval(.plot_and_jitter())
          
          # Update sample matrix and energy state
          # new_sm[wp, ] <- cbind(z = 1, covars[new_conf[wp, 1], ])# finite
          new_sm[wp, ] <- c(1, new_conf[wp, 2:3], covars[new_conf[wp, 1], all.vars(eqn)[-1]])
          new_energy <- data.frame(obj = .objMKV(
            eqn = eqn, sm = new_sm, covars = covars, vgm = vgm, krige.stat = krige.stat, debug.level = 0, 
            k = k, ...))
          
          # Avoid 'LDLfactor' error in 'krige' function
          # https://stat.ethz.ch/pipermail/r-sig-geo/2009-November/006919.html
          if (is.na(new_energy[1])) {
            new_energy <- old_energy
            new_conf <- old_conf
            new_sm <- old_sm
            message("\nskipped 'singular matrix' error in 'krige'-function")
          }

          # Evaluate the new system configuration
          accept <- .acceptSPSANN(old_energy[[1]], new_energy[[1]], actual_temp)
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
          if (track) energies[k, ] <- new_energy
          
          # Record best energy state
          if (new_energy[[1]] < best_energy[[1]] / 1.0000001) {
            best_k <- k
            best_conf <- new_conf
            best_energy <- new_energy
            best_old_energy <- old_energy
            old_conf <- old_conf
            best_sm <- new_sm
            best_old_sm <- old_sm
          }
          
          # Update progress bar
          eval(.update_progress())

        } # End loop through points

      } # End the chain

      # Check the proportion of accepted jitters in the first chain
      eval(.check_first_chain())
      
      # Count the number of chains without any change in the objective function.
      # Restart with the previously best configuration if it exists.
      if (n_accept == 0) {
        no_change <- no_change + 1
        if (no_change > schedule$stopping) {
          # if (new_energy > best_energy * 1.000001) {
            # old_conf <- old_conf
            # new_conf <- best_conf
            # old_energy <- best_old_energy
            # new_energy <- best_energy
            # new_sm <- best_sm
            # old_sm <- best_old_sm
            # no_change <- 0
            # cat("\nrestarting with previously best configuration\n")
          # } else { 
            break 
          # }
        }
        if (verbose) {
          cat("\n", no_change, "chain(s) with no improvement... stops at", schedule$stopping, "\n")
        }
      } else {
        no_change <-  0
      }
      
      # Update control parameters
      # Testing new parametes 'x_min0' and 'y_min0' (used with finite 'candi')
      actual_temp <- actual_temp * schedule$temperature.decrease
      x.max <- x_max0 - (i / schedule$chains) * (x_max0 - x.min) + cellsize[1] + x_min0
      y.max <- y_max0 - (i / schedule$chains) * (y_max0 - y.min) + cellsize[2] + y_min0
      
    } # End the annealing schedule
    
    # Prepare output
    eval(.prepare_output())
  }
# INTERNAL FUNCTION - CALCULATE THE ENERGY STATE VALUE ########################################################
# eqn: equation
# sm: sample matrix
# covars: covariates
# vgm: variogram model
# krige.stat: statistic to summarize the kriging variance
.objMKV <-
  function (eqn, sm, covars, vgm, krige.stat, k, ...) {
    
    # NOT ANYMORE!!!
    # We use 'set = list(cn_max = 1e10)' to avoid the LDFfactor error,
    # but do not accept the new system configuration.
    # https://stat.ethz.ch/pipermail/r-sig-geo/2009-November/006919.html
    # res <- gstat::krige(formula = eqn, locations = ~ x + y, data = sm,
                        # newdata = covars, model = vgm,
                        # set = list(cn_max = 1e10),
                        # ...)$var1.var
    # Error in predict.gstat(g, newdata = newdata, block = block, nsim = nsim, 
    # : gstat: value not allowed for: QRfactor not yet implemented 
    # I do not know the reason for this error, but it comes from using 
    # 'set = list(cn_max = 1e10)' above. I try to solve with 'tryCatch'!
    res <- NA
    try(
      res <- gstat::krige(
        formula = eqn, locations = ~ x + y, data = sm, newdata = covars, model = vgm, ...)$var1.var,
      silent = TRUE)
    
    # Calculate the energy state value
    if (krige.stat == "mean") { # Mean kriging variance
      res <- ifelse(k == 0, mean(res, na.rm = TRUE), mean(res))
    } else { # Maximum kriging variance
      res <- ifelse(k == 0, max(res, na.rm = TRUE), max(res))
    }
    
    # Output
    return (res)
  }
# INTERNAL FUNCTION - PREPARE STARTING SAMPLE MATRIX ##########################################################
# n_pts: number of points
# eqn: equation
# pts: points
# covars: covariates
.smMKV <-
  function (n_pts, eqn, pts, covars) {
    
    z <- rep(1, n_pts)
    if (stats::terms(eqn)[[3]] == 1) { # Simple and ordinary kriging
      # sm <- data.frame(z, pts[, 2:3])
      sm <- data.frame(z, pts[, c("x", "y")])
      colnames(sm) <- c("z", "x", "y")
      
    } else { # Universal kriging
      # sm <- data.frame(z, pts[, 2:3], covars[pts[, 1], all.vars(eqn)[-1]])
      sm <- data.frame(z, pts[, c("x", "y")], covars[pts[, 1], all.vars(eqn)[-1]])
      colnames(sm) <- c("z", "x", "y", all.vars(eqn)[-1])
    }
    
    # Output
    return (sm)
  }
# INTERNAL FUNCTION - PREPARE THE PREDICITON GRID #############################################################
# eqn: equation
# candi: prediction (evaluation) locations
# covars: covariates
.covarsMKV <-
  function (eqn, pred.grid, covars) {
    
    if (stats::terms(eqn)[[3]] == 1) { # Simple and ordinary kriging
      covars <- data.frame(pred.grid)
      colnames(covars) <- c("x", "y")
      
    } else { # Universal kriging
      covars <- as.data.frame(covars)
      covars <- data.frame(pred.grid, covars[, all.vars(eqn)[-1]])
      colnames(covars) <- c("x", "y", all.vars(eqn)[-1])
    }
    
    # Output
    return (covars)
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS #########################################################################
.checkMKV <-
  function (covars, eqn, vgm, krige.stat, candi, eval.grid) {
    
    # 'covars' compared to 'candi' e 'eval.grid'
    if (!missing(covars)) {
      covars_rows <- ifelse(is.vector(covars), length(covars), nrow(covars))
      if (!missing(eval.grid)) {
        if (covars_rows != nrow(eval.grid)) {
          res <- "'covars' and 'eval.grid' must have the same number of rows"
          return (res)
        }
      } else {
        if (covars_rows != nrow(candi)) {
          res <- "'covars' and 'candi' must have the same number of rows"
          return (res)
        }
      }
      # if (is.vector(covars)) {
      #   if (nrow(candi) != length(covars)) {
      #     res <- "'candi' and 'covars' must have the same number of rows"
      #     return (res)
      #   }
      # } else {
      #   if (nrow(candi) != nrow(covars)) {
      #     res <- "'candi' and 'covars' must have the same number of rows"
      #     return (res)
      #   }
      # }
    }
    
    # eqn
    if (missing(eqn)) {
      res <- "'eqn' must be a formula with dependent variable 'z'"
      return (res)
    } else {
      eqn_not_formula <- !inherits(eqn, "formula")
      eqn_dependent_not_z <- all.vars(eqn)[1] != "z"
      if (eqn_not_formula || eqn_dependent_not_z) {
        res <- "'eqn' must be a formula with dependent variable 'z'"
        return (res)
      }
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
      res <- paste("'krige.stat = ", krige.stat, "' is not supported", sep = "")
      return (res)
    }
  }
# FUNCTION - CANCLULATE THE OBJECTIVE FUNCTION VALUE ##########################################################
#' @export
#' @rdname optimMKV
objMKV <-
  function (points, candi, 
            # eval.grid,
            # MKV
            covars, eqn, vgm, krige.stat = "mean", ...) {
    
    # Check suggests
    pkg <- c("gstat")
    eval(.check_suggests())
    
    # Check other arguments
    check <- .checkMKV(
      covars = covars, eqn = eqn, vgm = vgm, krige.stat = krige.stat,
      # eval.grid = eval.grid, 
      candi = candi)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare prediction grid with covars
    # if (!missing(eval.grid)) { # Use coarser prediction (evaluation) grid
      # covars <- .covarsMKV(eqn = eqn, pred.grid = eval.grid, covars = covars)
    # } else { # Use candi as prediction (evaluation) grid
      covars <- .covarsMKV(eqn = eqn, pred.grid = candi[, c(2:3)], covars = covars)
    # }
    sm <- .smMKV(n_pts = n_pts, eqn = eqn, pts = points, covars = covars)
    
    # Energy state
    res <- .objMKV(eqn = eqn, sm = sm, covars = covars, vgm = vgm, krige.stat = krige.stat, k = 0, ...)
    
    # Output
    return (res)
  }
