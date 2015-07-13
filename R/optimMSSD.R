#' Optimization of sample configurations for spatial interpolation
#'
#' Optimize a sample configuration for spatial interpolation. The criterion 
#' used is the mean squared shortest distance (\bold{MSSD}).
#'
#' @template spJitter_doc
#' @template spSANN_doc
#' @template MOOP_doc
#'
#' @return
#' \code{optimMSSD} returns a matrix: the optimized sample configuration.
#'
#' \code{objMSSD} returns a numeric value: the energy state of the sample
#' configuration - the objective function value.
#' 
#' @note
#' This function was derived with modifications from the method known as the 
#' \emph{spatial coverage sampling} originally proposed by Brus, de Gruijter and
#' van Groenigen (2006), and implemented in the R-package 
#' \pkg{\link[spcosa]{spcosa}} by Dennis Walvoort, Dick Brus and Jaap de 
#' Gruijter.
#' 
#' @references
#' Brus, D. J.; de Gruijter, J. J.; van Groenigen, J. W. Designing spatial
#' coverage samples using the k-means clustering algorithm. In: P. Lagacherie,
#' A. M.; Voltz, M. (Eds.) \emph{Digital soil mapping - an introductory
#' perspective}. Elsevier, v. 31, p. 183-192, 2006.
#'
#' de Gruijter, J. J.; Brus, D.; Bierkens, M.; Knotters, M. \emph{Sampling for
#' natural resource monitoring}. Berlin: Springer, p. 332, 2006.
#'
#' Walvoort, D. J. J.; Brus, D. J.; de Gruijter, J. J. An R package for spatial
#' coverage sampling and random sampling from compact geographical strata by
#' k-means. \emph{Computers and Geosciences}. v. 36, p. 1261-1267, 2010.
#'
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso
#' \code{\link[raster]{distanceFromPoints}}, \code{\link[spcosa]{stratify}}.
#' @aliases optimMSSD objMSSD
#' @concept spatial interpolation
#' @export
#' @examples
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' set.seed(2001)
#' res <- optimMSSD(points = 100, candi = candi)
#' tail(attr(res, "energy.state"), 1) # 11531.03
#' objMSSD(candi = candi, points = res)
# FUNCTION - MAIN ##############################################################
optimMSSD <-
  function (points, candi, iterations = 100, 
    # SPSANN
    x.max, x.min, y.max, y.min,
    acceptance = list(initial = 0.99, cooling = iterations / 10),
    stopping = list(max.count = iterations / 10), plotit = FALSE, track = FALSE,
    boundary, progress = TRUE, verbose = FALSE, greedy = FALSE,
    # MOOP
    weights = NULL, nadir = NULL, utopia = NULL) {
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Calculate the initial energy state. The distance matrix is calculated
    # using the SpatialTools::dist2(). The function .objMSSD() does the
    # squaring internaly.
    dm <- SpatialTools::dist2(candi[, 2:3], conf0[, 2:3])
    energy0 <- .objMSSD(x = dm)
    
    # other settings for the simulated annealing algorithm
    old_dm <- dm
    best_dm <- dm
    count <- 0
    old_energy <- energy0
    best_energy <- Inf
    if (progress) {
      pb <- utils::txtProgressBar(min = 1, max = iterations, style = 3) 
    }
    time0 <- proc.time()
    
    # Begin the main loop
    for (k in 1:iterations) {
      
      # Plotting and jittering
      eval(.plot_and_jitter())
      
      # Update the matrix of distances and calculate the new energy state (Cpp)
      x2 <- matrix(new_conf[wp, 2:3], nrow = 1)
      new_dm <- .updateMSSDCpp(x1 = candi[, 2:3], x2 = x2, dm = old_dm, 
                               idx = wp)
      new_energy <- .objMSSD(new_dm)
      
      # Update the matrix of distances and calculate the new energy state (R)
      # x2 <- SpatialTools::dist2(coords = candi[, 2:3], coords2 = x2)
      # new_dm <- old_dm
      # new_dm[, wp] <- x2
      # new_energy <- mean(apply(new_dm, 1, min) ^ 2)
      
      # Evaluate the new system configuration
      if (greedy) {
        random_prob <- 1
      } else {
        random_prob <- stats::runif(1)
      }
      actual_prob <- acceptance[[1]] * exp(-k / acceptance[[2]])
      if (track) accept_probs[k] <- actual_prob
      if (new_energy <= old_energy) {
        old_conf <- new_conf
        old_energy <- new_energy
        count <- 0
        old_dm <- new_dm
      } else {
        if (new_energy > old_energy & random_prob <= actual_prob) {
          old_conf <- new_conf
          old_energy <- new_energy
          count <- count + 1
          old_dm <- new_dm
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... p = ",
                random_prob, "\n")
          }
          } else {
            new_energy <- old_energy
            new_conf <- old_conf
            count <- count + 1
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
        best_dm <- new_dm
        best_old_dm <- old_dm
      }

      # Freezing parameters
      if (count == stopping[[1]]) {
        if (new_energy > best_energy * 1.000001) {
          old_conf <- old_conf
          new_conf <- best_conf
          old_energy <- best_old_energy
          new_energy <- best_energy
          count <- 0
          new_dm <- best_dm
          old_dm <- best_old_dm
          cat("\n", "reached maximum count with suboptimal configuration\n")
          cat("\n", "restarting with previously best configuration\n")
          cat("\n", count, "iteration(s) with no improvement... stops at",
              stopping[[1]], "\n")
        } else {
          break
        }
      }
      if (progress) utils::setTxtProgressBar(pb, k)
    }
    
    # Prepare output
    eval(.prepare_output())
  }
# FUNCTION - CALCULATE THE CRITERION VALUE #####################################
#' @rdname optimMSSD
#' @export
objMSSD <-
  function (points, candi) {
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Compute the matrix of distances
    dm <- SpatialTools::dist2(coords1 = candi[, 2:3], coords2 = points[, 2:3])
    
    # Calculate the energy state
    res <- .objMSSD(dm)
    
    # Output
    return (res)
  }
# INTERNAL FUNCTION - NEAREST POINT ############################################
# Function to identify the nearest point.
.nearestPoint <-
  function (dist.mat, which.pts) {
    if (dim(dist.mat)[1] != dim(dist.mat)[2] || !is.matrix(dist.mat)) {
      stop ("'dist.mat' should be a n x n matrix")
    }
    if (!is.vector(which.pts) || length(which.pts) >= dim(dist.mat)[2]) {
      stop (paste("'which.pts' should be a vector of length smaller than ",
                  dim(dist.mat)[2], sep = ""))
    }
    sub_dm <- dist.mat[, which.pts]
    min_dist_id <- apply(sub_dm, MARGIN = 1, FUN = which.min)
    return (min_dist_id)
  }
# End!
