#' Optimization of sample configurations for spatial interpolation
#'
#' Optimize a sample configuration for spatial interpolation. The criterion 
#' used is
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
#' \code{optimMSSD} returns a matrix: the optimized sample configuration with
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
#' require(pedometrics)
#' require(sp)
#' require(SpatialTools)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' set.seed(2001)
#' res <- optimMSSD(points = 100, candi = candi, iterations = 100)
#' tail(attr(res, "energy.state"), 1) # 11531.03
#' objMSSD(candi = candi, points = res)
# FUNCTION - MAIN ##############################################################
optimMSSD <-
  function (points, candi, x.max, x.min, y.max, y.min, iterations = 10000,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            boundary, progress = TRUE, verbose = TRUE, greedy = FALSE,
            weights, nadir, utopia) {
    
    # Check spsann arguments ###################################################
    check_spsann_arguments <- 
      function (...) {parse(text = readLines("tools/check-spsann-arguments.R"))}
    eval(check_spsann_arguments())
    ############################################################################
    
    # Set plotting options ####################################################
    plotting_options <- 
      function (...) {parse(text = readLines("tools/plotting-options.R"))}
    eval(plotting_options())
    ############################################################################
    
    # Prepare sample points
    #n_candi <- nrow(candi)
    #points <- .spsannPoints(points = points, candi = candi, n.candi = n_candi)
    #n_pts <- nrow(points)
    #conf0 <- points
    #old_conf <- conf0
    # Prepare points and candi #################################################
    prepare_points <- 
      function (...) {parse(text = readLines("tools/prepare-points.R"))}
    eval(prepare_points())
    ############################################################################
    
    # Prepare for jittering ####################################################
    prepare_jittering <- 
      function (...) {parse(text = readLines("tools/prepare-jittering.R"))}
    eval(prepare_jittering())
    ############################################################################
    
    # Calculate the initial energy state. The distance matrix is calculated
    # using the SpatialTools::dist2(). The function .calcMSSDCpp() does the
    # squaring internaly.
    # ASR: write own distance function in C++
    dm <- SpatialTools::dist2(candi[, 2:3], conf0[, 2:3])
    energy0 <- .calcMSSDCpp(x = dm)
    #energy0 <- mean(apply(dm, 1, min) ^ 2)
    
    # other settings for the simulated annealing algorithm
    MOOP = FALSE
    old_dm        <- dm
    best_dm       <- dm
    count         <- 0
    old_energy    <- energy0
    best_energy   <- Inf
    energy_states <- vector()
    accept_probs  <- vector()
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()
    
    # Begin the main loop
    for (k in 1:iterations) {
      
      # Jitter one of the points and update x.max and y.max
      # ASR: spJitterFinite() can be improved implementing it in C++
      wp <- sample(c(1:n_pts), 1)
      new_conf <- spJitterFinite(old_conf, candi, x.max, x.min, y.max,
                                 y.min, wp)
      x.max <- x_max0 - (k / iterations) * (x_max0 - x.min)
      y.max <- y_max0 - (k / iterations) * (y_max0 - y.min)
      
      # Update the distance matrix and calculate the new energy state
      x2 <- matrix(new_conf[wp, 2:3], nrow = 1)
      
      # Update the matrix of distances in C++
      new_dm <- .updateMSSDCpp(x1 = candi[, 2:3], x2 = x2, dm = old_dm, 
                               idx = wp)
      
      # Update the matrix of distances in R
      #x2 <- SpatialTools::dist2(coords = candi[, 2:3], coords2 = x2)
      #new_dm <- old_dm
      #new_dm[, wp] <- x2
      
      # new_energy <- mean(apply(new_dm, 1, min) ^ 2)
      new_energy <- .calcMSSDCpp(new_dm)
      
      # ASR: This is to test the 'update' function
      #a <- objMSSD(candi = candi, points = new_conf)
      #if (round(new_energy, 2) != round(a, 2)) {
      # print(new_energy)
      # print(a)
      # break
      #}
            
      # Evaluate the new system configuration
      if (greedy) {
        random_prob <- 1
      } else {
        random_prob <- runif(1)
      }
      actual_prob <- acceptance[[1]] * exp(-k / acceptance[[2]])
      accept_probs[k] <- actual_prob
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
            new_conf   <- old_conf
            count      <- count + 1
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
        best_conf     <- new_conf
        best_energy     <- new_energy
        best_old_energy <- old_energy
        old_conf      <- old_conf
        best_dm         <- new_dm
        best_old_dm     <- old_dm
      }

      # Plotting
      #if (any(round(seq(1, iterations, 10)) == k)) {
      if (plotit && pedometrics::is.numint(k / 10)) {
        if (plotit){
          .spSANNplot(energy0, energy_states, k, acceptance,
                      accept_probs, boundary, new_conf[, 2:3],
                      conf0[, 2:3], y_max0, y.max, x_max0, x.max,
                      best.energy = best_energy, best.k = best_k, MOOP = MOOP)
        }
      }

      # Freezing parameters
      if (count == stopping[[1]]) {
        if (new_energy > best_energy * 1.000001) {
          old_conf <- old_conf
          new_conf <- best_conf
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
    res <- .spSANNout(new_conf, energy0, energy_states, time0, MOOP = MOOP)
    return (res)
  }
# FUNCTION - CALCULATE THE CRITERION VALUE #####################################
#' @rdname optimMSSD
#' @export
objMSSD <-
  function (candi, points) {
    
    # Prepare points and candi #################################################
    prepare_points <- 
      function (...) {parse(text = readLines("tools/prepare-points.R"))}
    eval(prepare_points())
    ############################################################################
    
    coords1 <- candi[, 2:3]
    coords2 <- points[, 2:3]
    dm <- SpatialTools::dist2(coords1 = coords1, coords2 = coords2)
    res <- mean(apply(dm, 1, min) ^ 2)
    #res <- .calcMSSDCpp(dm)
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
