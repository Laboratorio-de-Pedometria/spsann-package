#' Optimization of sample configurations using a user-defined objective function
#'
#' Optimize a sample configuration using a user-defined objective function.
#' 
#' @template spJitter_doc
#' @template spSANN_doc
#' 
#' @param fun A function defining the objective function that should be
#' used to evaluate the energy state of the system configuration at each 
#' iteration. See \sQuote{Details} for more information.
#' 
#' @param ... Other arguments passed to the objective function. See 
#' \sQuote{Details} for more information.
#'
#' @details
#' The user-defined objective function has to be an object of class 
#' \link[base]{function}. It has to include the argument \code{points}, which is
#' defined internally as a matrix with three columns: \code{[, 1]} the 
#' identification of each sample point (1, 2, ..., n), \code{[, 2]} the 
#' x-coordinates, and \code{[, 3]} the y-coordinates. The identification is 
#' useful to retrieve information from any data matrix used by the objective
#' function defined by the user.
#' 
#' @return
#' \code{optimUSER} returns a matrix: the optimized sample configuration with
#' the evolution of the energy state during the optimization as an attribute.
#'
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @aliases optimUSER
#' @keywords spatial optimize
#' @concept simulated annealing
#' @export
#' @examples
#' require(pedometrics)
#' require(sp)
#' require(rgeos)
#' require(SpatialTools)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' coordinates(candi) <- ~ x + y
#' gridded(candi) <- TRUE
#' boundary <- as(candi, "SpatialPolygons")
#' boundary <- gUnionCascaded(boundary)
#' candi <- coordinates(candi)
#' candi <- matrix(cbind(1:nrow(candi), candi), ncol = 3)
#' x.max <- diff(bbox(boundary)[1, ])
#' y.max <- diff(bbox(boundary)[2, ])
#' 
#' # Define the objective function - number of points per lag distance class
#' objUSER <-
#'   function (points, lags, n_lags, n_pts) {
#'     dm <- SpatialTools::dist1(points[, 2:3])
#'     ppl <- vector()
#'     for (i in 1:n_lags) {
#'       n <- which(dm > lags[i] & dm <= lags[i + 1], arr.ind = TRUE)
#'       ppl[i] <- length(unique(c(n)))
#'     }
#'     distri <- rep(n_pts, n_lags)
#'     res <- sum(distri - ppl)
#'   }
#' lags <- seq(1, 1000, length.out = 10)
#' 
#' # Run the optimization using the user-defined objective function
#' set.seed(2001)
#' timeUSER <- Sys.time()
#' resUSER <- optimUSER(points = 100, fun = objUSER, lags = lags, 
#'                      n_lags = 9, n_pts = 100,
#'                      candi = candi, x.max = x.max, x.min = 40, y.max = y.max,
#'                      y.min = 40, boundary = boundary, iterations = 1000)
#' timeUSER <- Sys.time() - timeUSER
#' 
#' # Run the optimization using the respective function implemented in spsann
#' set.seed(2001)
#' timePPL <- Sys.time()
#' resPPL <- optimPPL(points = 100, candi = candi, lags = lags,  
#'                    criterion = "distribution", x.max = x.max, x.min = 40, 
#'                    y.max = y.max, y.min = 40, boundary = boundary,
#'                    iterations = 1000)
#' timePPL <- Sys.time() - timePPL
#' 
#' # Compare results
#' timeUSER
#' timePPL
#' lapply(list(resUSER, resPPL), countPPL, lags = lags, pairs = FALSE)
#' x <- attr(resUSER, "energy.state")
#' y <- attr(resPPL, "energy.state")
#' sapply(list(x, y), tail, 1)
#' plot(x, y, asp = 1)
#' abline(0, 1, col = "red")
# FUNCTION - MAIN ##############################################################
optimUSER <-
  function (points, fun, ...,
            candi, x.max, x.min, y.max, y.min, iterations = 10000,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            boundary, progress = TRUE, verbose = TRUE, greedy = FALSE) {
    
    # Check arguments
    # http://www.r-bloggers.com/a-warning-about-warning/
    check <- .spSANNcheck(points = points, candi = candi, x.max = x.max, 
                          x.min = x.min, y.max = y.max, y.min = y.min, 
                          iterations = iterations, acceptance = acceptance,
                          stopping = stopping, plotit = plotit, 
                          boundary = boundary, progress = progress, 
                          verbose = verbose)
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
    
    # Initial energy state
    energy0 <- .energyState(fun = fun, points = old_conf, ...)
    
    # Other settings for the simulated annealing algorithm
    MOOP <- FALSE
    count <- 0
    old_energy <- energy0
    best_energy <- Inf
    energies <- vector()
    accept_probs <- vector()
    x_max0 <- x.max
    y_max0 <- y.max
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()
    
    # Begin the iterations    
    for (k in 1:iterations) {
      
      # Jitter one of the points and update x.max and y.max
      wp <- sample(1:n_pts, 1)
      new_conf <- spJitterFinite(points = old_conf, candi = candi,
                                 x.max = x.max, x.min = x.min, y.max = y.max,
                                 y.min = y.min, which.point = wp)
      x.max <- x_max0 - (k / iterations) * (x_max0 - x.min)
      y.max <- y_max0 - (k / iterations) * (y_max0 - y.min)
      
      # New energy state
      new_energy <- .energyState(fun = fun, points = new_conf, ...)
      
      # Evaluate the new system configuration
      if (greedy) {
        random_prob <- 1
      } else {
        random_prob <- runif(1)
      }
      actual_prob <- acceptance$initial * exp(-k / acceptance$cooling)
      accept_probs[k] <- actual_prob
      if (new_energy <= old_energy) {
        # Always accepts a better energy
        old_conf <- new_conf
        old_energy <- new_energy
        count <- 0
      } else {
        if (new_energy > old_energy && random_prob <= actual_prob) {
          # Accepts a worse energy depending on the probability
          old_conf <- new_conf
          old_energy <- new_energy
          count <- count + 1
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
                stopping$max.count, "\n")
          }
        }
      }
      
      # Best energy state
      energies[k] <- new_energy
      if (new_energy < best_energy / 1.0000001) {
        best_k <- k
        best_conf <- new_conf
        best_energy <- new_energy
        best_old_energy <- old_energy
        old_conf <- old_conf
      }
      
      # Plotting
      if (plotit && pedometrics::is.numint(k / 10)) {
        .spSANNplot(energy0 = energy0, energies = energies, k = k, 
                    acceptance = acceptance, accept_probs = accept_probs,
                    boundary = boundary, new_conf = new_conf[, 2:3], 
                    conf0 = conf0[, 2:3], y_max0 = y_max0, y.max = y.max,
                    x_max0 = x_max0, x.max = x.max, best.energy = best_energy,
                    best.k = best_k, MOOP = MOOP)
      }
      
      # Freezing parameters
      if (count == stopping$max.count) {
        if (new_energy > best_energy * 1.000001) {
          old_conf <- old_conf
          new_conf <- best_conf
          old_energy <- best_old_energy
          new_energy <- best_energy
          count <- 0
          energies[k] <- new_energy
          cat("\n", "reached maximum count with suboptimal configuration\n")
          cat("\n", "restarting with previously best configuration\n")
          cat("\n", count, "iteration(s) with no improvement... stops at",
              stopping$max.count, "\n")
        } else {
          break
        }
      }
      if (progress) setTxtProgressBar(pb, k)
    }
    if (progress) close(pb)
    res <- .spSANNout(new_conf = new_conf, energy0 = energy0, MOOP = MOOP,
                      energies = energies, time0 = time0)
    return (res)
  }
# INTERNAL FUNCTION - CALCULATE DE ENERGY STATE ################################
.energyState <- 
  function (fun, points, ...) {
    if (missing(fun) || missing(points)) {
      stop ("'fun' and 'points' are mandatory arguments")
    }
    return (do.call(fun, list(points, ...)))
  }
