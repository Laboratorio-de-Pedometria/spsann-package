#' Optimization of sample configurations for spatial interpolation (I)
#'
#' Optimize a sample configuration for spatial interpolation with a 'known' auto- or cross-correlation model, 
#' e.g. simple and ordinary (co)kriging. The criterion used is the mean squared shortest distance (__MSSD__) 
#' between sample locations and prediction locations.
#'
# @inheritParams spJitter
#' @template spSANN_doc
#' @template spJitter_doc
#' @template schedule_doc
#' 
#' @param eval.grid (Experimental) Data frame or matrix with the objective function evaluation locations. Like
#' `candi`, `eval.grid` must have two columns in the following order: `[, "x"]`, the projected x-coordinates,
#' and `[, "y"]`, the projected y-coordinates.
#' 
#' @details 
#' \subsection{Mean squared shortest distance}{
#' This objective function is based on the knowledge that the simple and ordinary (co)kriging prediction error 
#' variance only depends upon the separation distance between sample locations: the larger the distance, the
#' larger the prediction error variance. As such, the better the spread of the sample locations in the spatial
#' domain, the smaller the overall simple/ordinary (co)kriging prediction error variance. This is the purpose 
#' of using a regular grid of sample locations.
#' 
#' However, a regular grid usually is suboptimal, especially if the spatial domain is irregularly shaped. Thus
#' the need for optimization, that is based on measuring the goodness of the spread of sample locations in the
#' spatial domain. To measure this spread we can compute the distance from every sample location to each of the
#' prediction locations placed on a fine grid covering the entire spatial domain. Next, for every prediction
#' location we find the closest sample location and record its distance. The mean of these squared distances
#' over all prediction location will measure the spread of the sample locations.
#' 
#' During the optimization, we try to reduce this measure -- the mean squared shortest distance -- between
#' sample and prediction locations. (This is also know as _spatial coverage sampling_, see the R-package
#' __[spcosa](https://CRAN.R-project.org/package=spcosa)__.)
#' }
#'
#' @return
#' \code{optimMSSD} returns an object of class \code{OptimizedSampleConfiguration}: the optimized sample
#' configuration with details about the optimization.
#'
#' \code{objMSSD} returns a numeric value: the energy state of the sample configuration -- the objective
#' function value in square map units, generally m^2 or km^2.
#' 
#' @note
#' \subsection{Sample configuration for spatial interpolation}{
#' A sample configuration optimized for spatial interpolation such as simple and ordinary (co)kriging is not
#' necessarily appropriate for estimating the parameters of the spatial autocorrelation model, i.e. the
#' parameters of the variogram model. See \code{\link[spsann]{optimPPL}} for more information on the 
#' optimization of sample configurations for variogram identification and estimation.
#' }
#' 
#' @references
#' Brus, D. J.; de Gruijter, J. J.; van Groenigen, J.-W. Designing spatial coverage samples using the k-means
#' clustering algorithm. In: P. Lagacherie,A. M.; Voltz, M. (Eds.) _Digital soil mapping -- an introductory
#' perspective_. Elsevier, v. 31, p. 183-192, 2006.
#'
#' de Gruijter, J. J.; Brus, D.; Bierkens, M.; Knotters, M. _Sampling for natural resource monitoring_.
#' Berlin: Springer, p. 332, 2006.
#'
#' Walvoort, D. J. J.; Brus, D. J.; de Gruijter, J. J. An R package for spatial coverage sampling and random 
#' sampling from compact geographical strata by k-means. _Computers and Geosciences_. v. 36, p. 
#' 1261-1267, 2010.
#'
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso
#' \code{[distanceFromPoints](https://CRAN.R-project.org/package=raster)},
#' \code{[stratify](https://CRAN.R-project.org/package=spcosa)}.
#' @aliases optimMSSD objMSSD MSSD
#' @concept spatial interpolation
#' @export
#' @examples
#' #####################################################################
#' # NOTE: The settings below are unlikely to meet your needs.         #
#' #####################################################################
#' data(meuse.grid, package = 'sp')
#' candi <- meuse.grid[, 1:2]
#' schedule <- scheduleSPSANN(
#'   chains = 1, initial.temperature = 5000000,
#'   x.max = 1540, y.max = 2060, x.min = 0, y.min = 0, cellsize = 40)
#' set.seed(2001)
#' res <- optimMSSD(points = 10, candi = candi, schedule = schedule)
#' data.frame(
#'   expected = 247204.8,
#'   objSPSANN = objSPSANN(res),
#'   objMSSD = objMSSD(candi = candi, points = res)
#' )
#' 
# FUNCTION - MAIN #############################################################################################
optimMSSD <-
  function (points, candi, eval.grid,
            # SPSANN
            schedule = scheduleSPSANN(), plotit = FALSE, track = FALSE, 
            boundary, progress = "txt", verbose = FALSE) {
    
    # Objective function name
    objective <- "MSSD"
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Calculate the initial energy state. The distance matrix is calculated using the SpatialTools::dist2(). 
    # The function .objMSSD() does the squaring internaly.
    # 'old_conf' is used instead of 'conf0' because the former holds information on both fixed and free points.
    # dm <- SpatialTools::dist2(candi[, 2:3], conf0[, 2:3])
    if (!missing(eval.grid)) {
      dm <- SpatialTools::dist2(coords1 = eval.grid, coords2 = old_conf[, 2:3])
    } else {
      dm <- SpatialTools::dist2(coords1 = candi[, 2:3], coords2 = old_conf[, 2:3])
    }
    energy0 <- data.frame(obj = .objMSSD(x = dm))
    
    # Other settings for the simulated annealing algorithm
    old_dm <- dm
    best_dm <- dm
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
          
          # Update the matrix of distances and calculate the new energy state
          x2 <- matrix(new_conf[wp, 2:3], nrow = 1)
          if (!missing(eval.grid)) {
            new_dm <- .updateMSSDCpp(x1 = eval.grid, x2 = x2, dm = old_dm, idx = wp)
          } else {
            new_dm <- .updateMSSDCpp(x1 = candi[, 2:3], x2 = x2, dm = old_dm, idx = wp)
          }
          new_energy <- data.frame(obj = .objMSSD(new_dm))
          
          # Update the matrix of distances and calculate the new energy state
          # x2 <- SpatialTools::dist2(coords = candi[, 2:3], coords2 = x2)
          # new_dm <- old_dm
          # new_dm[, wp] <- x2
          # new_energy <- mean(apply(new_dm, 1, min) ^ 2)
          
          # Evaluate the new system configuration
          accept <- .acceptSPSANN(old_energy[[1]], new_energy[[1]], actual_temp)
          if (accept) {
            old_conf <- new_conf
            old_energy <- new_energy
            old_dm <- new_dm
            n_accept <- n_accept + 1
          } else {
            new_energy <- old_energy
            new_conf <- old_conf
            # new_dm <- old_dm
          }
          if (track) energies[k, ] <- new_energy
          
          # Record best energy state
          if (new_energy[[1]] < best_energy[[1]] / 1.0000001) {
            best_k <- k
            best_conf <- new_conf
            best_energy <- new_energy
            best_old_energy <- old_energy
            old_conf <- old_conf
            best_dm <- new_dm
            best_old_dm <- old_dm
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
          # if (new_energy[[1]] > best_energy[[1]] * 1.000001) {
            # old_conf <- old_conf
            # new_conf <- best_conf
            # old_energy <- best_old_energy
            # new_energy <- best_energy
            # new_dm <- best_dm
            # old_dm <- best_old_dm
            # no_change <- 0
            # new_sm <- best_sm
            # old_sm <- best_old_sm
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
      actual_temp <- actual_temp * schedule$temperature.decrease
      x.max <- x_max0 - (i / schedule$chains) * (x_max0 - x.min) + cellsize[1]
      y.max <- y_max0 - (i / schedule$chains) * (y_max0 - y.min) + cellsize[2]
      
    } # End the annealing schedule
    
    # Prepare output
    eval(.prepare_output())
    return (res)
  }
# FUNCTION - CALCULATE THE CRITERION VALUE ####################################################################
#' @rdname optimMSSD
#' @export
objMSSD <-
  function (points, candi, eval.grid) {
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Compute the matrix of distances
    if (!missing(eval.grid)) {
      dm <- SpatialTools::dist2(coords1 = eval.grid, coords2 = points[, 2:3])
    } else {
      dm <- SpatialTools::dist2(coords1 = candi[, 2:3], coords2 = points[, 2:3])
    }
    
    # Calculate the energy state
    res <- .objMSSD(dm)
    
    # Output
    return (res)
  }
# INTERNAL FUNCTION - NEAREST POINT ###########################################################################
# Function to identify the nearest point.
.nearestPoint <-
  function (dist.mat, which.pts) {
    if (dim(dist.mat)[1] != dim(dist.mat)[2] || !is.matrix(dist.mat)) {
      stop ("'dist.mat' should be a n x n matrix")
    }
    if (!is.vector(which.pts) || length(which.pts) >= dim(dist.mat)[2]) {
      stop (paste("'which.pts' should be a vector of length smaller than ", dim(dist.mat)[2], sep = ""))
    }
    sub_dm <- dist.mat[, which.pts]
    min_dist_id <- apply(sub_dm, MARGIN = 1, FUN = which.min)
    return (min_dist_id)
  }
# INTERNAL FUNCTION ###########################################################################################
.subsetCandi <-
  function (candi, cellsize) {
    x <- seq(min(candi[, 1]), max(candi[, 1]), by = cellsize * 2)
    y <- seq(min(candi[, 2]), max(candi[, 2]), by = cellsize * 2)
    xy <- expand.grid(x, y)
    colnames(xy) <- colnames(candi)
    dup <- rbind(xy, candi)
    dup <- which(duplicated(dup)) - nrow(xy)
    return (dup)
  }
# End!
