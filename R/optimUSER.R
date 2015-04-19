#' Optimization of sample patterns using a user-defined objective function
#'
#' Optimize a sample pattern using a user-defined objective function.
#' 
#' @template spJitter_doc
#' @template spSANN_doc
#' 
#' @param
#'
#' @details
#' 
#' @return
#' \code{optimUSER} returns a matrix: the optimized sample pattern with
#' the evolution of the energy state during the optimization as an attribute.
#'
#' @references
#'
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @aliases optimUSER
#' @keywords spatial optimize
#' @concept simulated annealing
#' @export
#' @examples
#' 
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
    energy0 <- .energyState(fun = fun, points = old_conf[, 2:3], ...)
    
    # Other settings for the simulated annealing algorithm
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
      new_energy <- .energyState(fun = fun, points = new_conf[, 2:3], ...)
      
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
                    best.k = best_k)
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
    res <- .spSANNout(new_conf = new_conf, energy0 = energy0, 
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
