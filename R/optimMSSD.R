#' Optimization of sample patterns for spatial interpolation
#'
#' Optimize a sample pattern for spatial interpolation. The criterion used is
#' the mean squared shortest distance (\code{optimMSSD}). \code{objMSSD}
#' computes the MSSD between a set of points and all grid cells.
#'
#' @template spJitter_doc
#' @template spSANN_doc
#'
#' @details
#' Euclidean distances between points are calculated. This computation requires
#' the coordinates to be projected. The user is responsible for making sure that
#' this requirement is attained.
#' @return
#' \code{objMSSD} returns a numeric value: the MSSD between a set of points and
#' all grid cells.
#'
#' \code{optimMSSD} returns a matrix: the optimized sample pattern with
#' the evolution of the energy state during the optimization as an attribute.
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
#' @keywords spatial optimize
#' @aliases optimMSSD objMSSD
#' @concept simulated annealing
#' @importFrom fields rdist
#' @export
#' @examples
#' require(sp)
#' data(meuse.grid)
#' meuse.grid <- as.matrix(meuse.grid[, 1:2])
#' meuse.grid <- matrix(cbind(c(1:dim(meuse.grid)[1]), meuse.grid), ncol = 3)
#' points <- sample(c(1:dim(meuse.grid)[1]), 155)
#' points <- meuse.grid[points, ]
#' objMSSD(meuse.grid, points)
# FUNCTION - MAIN ##############################################################
optimMSSD <-
  function (points, candidates, x.max, x.min, y.max, y.min, iterations = 10000,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            boundary, progress = TRUE, verbose = TRUE) {
    
    check <- .spSANNcheck(points = points, candidates = candidates, 
                          x.max = x.max, x.min = x.min, y.max = y.max, 
                          y.min = y.min, iterations = iterations, 
                          acceptance = acceptance, stopping = stopping, 
                          plotit = plotit, boundary = boundary, 
                          progress = progress, verbose = verbose)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    if (plotit) {
      par0 <- par()
      on.exit(suppressWarnings(par(par0)))
    }
    n_candi <- nrow(candidates)
    if (is.integer(points) || pedometrics::is.numint(points)) {
      n_pts <- points
      points <- sample(1:n_candi, n_pts)
      points <- candidates[points, ]
    } else {
      n_pts <- nrow(points)
    }
    config0 <- points
    old_config <- config0

    # Calculate the initial energy state. The distance matrix is calculated
    # using the fields::rdist(). The function .calcMSSDCpp() does the squaring
    # internaly.
    # ASR: write own distance function in C++
    dm <- fields::rdist(candidates[, 2:3], config0[, 2:3])
    energy0 <- .calcMSSDCpp(dm)
    
    # other settings for the simulated annealing algorithm
    old_dm        <- dm
    new_dm        <- dm
    best_dm       <- dm
    count         <- 0
    old_energy    <- energy0
    best_energy   <- Inf
    energy_states <- vector()
    accept_probs  <- vector()
    x_max0        <- x.max
    y_max0        <- y.max
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()
    
    # Begin the main loop
    for (k in 1:iterations) {
      
      # Jitter one of the points and update x.max and y.max
      # ASR: spJitterFinite() can be improved implementing it in C++
      wp <- sample(c(1:n_pts), 1)
      new_config <- spJitterFinite(old_config, candidates, x.max,
                                       x.min, y.max, y.min, wp)
      x.max <- x_max0 - (k / iterations) * (x_max0 - x.min)
      y.max <- y_max0 - (k / iterations) * (y_max0 - y.min)
      
      # Update the distance matrix and calculate the new energy state
      x2 <- matrix(new_config[wp, 2:3], nrow = 1)
      new_dm <- .updateMSSDCpp(x1 = candidates[, 2:3], x2 = x2, dm = old_dm, 
                               idx = wp)
      new_energy <- .calcMSSDCpp(new_dm)
      
      # Evaluate the new system configuration
      random_prob     <- runif(1)
      actual_prob     <- acceptance[[1]] * exp(-k / acceptance[[2]])
      accept_probs[k] <- actual_prob
      if (new_energy <= old_energy) {
        old_config <- new_config
        old_energy <- new_energy
        count      <- 0
        old_dm     <- new_dm
      } else {
        if (new_energy > old_energy & random_prob <= actual_prob) {
          old_config <- new_config
          old_energy <- new_energy
          count      <- count + 1
          old_dm     <- new_dm
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... p = ",
                random_prob, "\n")
          }
          } else {
            new_energy <- old_energy
            new_config <- old_config
            count      <- count + 1
            new_dm     <- old_dm
            if (verbose) {
              cat("\n", count, "iteration(s) with no improvement... stops at",
                  stopping[[1]], "\n")
            }
          }
      }

      # Best energy state
      energy_states[k] <- new_energy
      if (new_energy < best_energy / 1.0000001) {
        best_k          <- k
        best_config     <- new_config
        best_energy     <- new_energy
        best_old_energy <- old_energy
        old_config      <- old_config
        best_dm         <- new_dm
        best_old_dm     <- old_dm
      }

      # Plotting
      if (any(round(seq(1, iterations, 10)) == k)) {
        if (plotit){
          .spSANNplot(energy0, energy_states, k, acceptance,
                      accept_probs, boundary, new_config[, 2:3],
                      config0[, 2:3], y_max0, y.max, x_max0, x.max)
        }
      }

      # Freezing parameters
      if (count == stopping[[1]]) {
        if (new_energy > best_energy * 1.000001) {
          old_config <- old_config
          new_config <- best_config
          old_energy <- best_old_energy
          new_energy <- best_energy
          count      <- 0
          new_dm     <- best_dm
          old_dm     <- best_old_dm
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
    if (plotit) par(par0)
    res <- .spSANNout(new_config, energy0, energy_states, time0)
    return (res)
  }
# FUNCTION - CALCULATE THE CRITERION VALUE #####################################
#' @rdname optimMSSD
#' @export
objMSSD <-
  function (candidates, points) {
    dm <- fields::rdist(candidates[, 2:3], points[, 2:3])
    res <- .calcMSSDCpp(dm)
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
