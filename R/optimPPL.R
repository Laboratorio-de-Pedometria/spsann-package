#' Optimization of sample configurations for variogram identification and 
#' estimation
#'
#' Optimize a sample configuration for variogram identification and estimation. 
#' A criterion is defined so that the optimized sample configuration has a 
#' given number of points or point-pairs contributing to each lag-distance 
#' class (\bold{PPL}).
#' 
#' @template spJitter_doc
#' @template spSANN_doc
#' @template MOOP_doc
#' @template PPL_doc
#' 
#' @return
#' \code{optimPPL} returns a matrix: the optimized sample configuration.
#'
#' \code{objPPL} returns a numeric value: the energy state of the sample 
#' configuration - the objective function value.
#'
#' \code{countPPL} returns a data.frame with three columns: a) the lower and b)
#' upper limits of each lag-distance class, and c) the number of points or 
#' point-pairs per lag-distance class.
#'
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @aliases optimPPL countPPL objPPL
#' @export
#' @examples
#' \dontrun{
#' # This example takes more than 5 seconds to run!
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' set.seed(2001)
#' res <- optimPPL(points = 100, candi = candi)
#' objSPSANN(res) # 160
#' objPPL(points = res, candi = candi)
#' countPPL(points = res, candi = candi)
#' }
# FUNCTION - MAIN ##############################################################
optimPPL <-
  function (points, candi, iterations = 100, 
    # PPL
    lags = 7, lags.type = "exponential", lags.base = 2, cutoff, 
    criterion = "distribution", distri, pairs = FALSE,
    # SPSANN
    x.max, x.min, y.max, y.min,
    acceptance = list(initial = 0.90, cooling = iterations / 10, 
                      by = "iterations", temperature = 5, calibrate = TRUE),
    stopping = list(max.count = iterations / 10), plotit = FALSE, track = FALSE,
    boundary, progress = TRUE, verbose = FALSE, greedy = FALSE,
    # MOOP
    weights = NULL, nadir = NULL, utopia = NULL) {
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Check other arguments
    check <- 
      .checkPPL(lags = lags, lags.type = lags.type, pairs = pairs,
                lags.base = lags.base, cutoff = cutoff, criterion = criterion,
                distri = distri, fun = "optimPPL")
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Prepare cutoff and lags
    cutoff <- .cutoffPPL(cutoff = cutoff, x.max = x.max, y.max = y.max)
    lags <- .lagsPPL(lags = lags, lags.type = lags.type, cutoff = cutoff, 
                     lags.base = lags.base)
    n_lags <- length(lags) - 1
    
    # Initial energy state: points or point-pairs
    dm <- SpatialTools::dist1(conf0[, 2:3])
    ppl <- .getPPL(lags = lags, n.lags = n_lags, dist.mat = dm, pairs = pairs)
    distri <- .distriPPL(n.lags = n_lags, n.pts = n_pts, criterion = criterion,
                         distri = distri, pairs = pairs)
    energy0 <- .objPPL(ppl = ppl, n.lags = n_lags, n.pts = n_pts,
                       criterion = criterion, distri = distri, pairs = pairs)
    
    # Other settings for the simulated annealing algorithm
    old_dm <- dm
    best_dm <- dm
    count <- 0
    old_energy <- energy0
    best_energy <- Inf
    if (progress) { 
      pb <- utils::txtProgressBar(min = 1, max = iterations, style = 3) 
    }
    time0 <- proc.time()
    
    # Begin the iterations
    k <- 0
    repeat {
      k <- k + 1
      
      # Plotting and jittering
      eval(.plot_and_jitter())
      
      # Update the distance matrix using a Cpp function
      new_dm <- .updatePPLCpp(x = new_conf[, 2:3], dm = old_dm, idx = wp)
      
      # Recalculate the full distance matrix
      #new_dm <- SpatialTools::dist1(coords = new_conf[, 2:3])
      
      # Update the distance matrix in R
      #x2 <- matrix(new_conf[wp, 2:3], nrow = 1)
      #x2 <- SpatialTools::dist2(coords = new_conf[, 2:3], coords2 = x2)
      #new_dm <- old_dm
      #new_dm[wp, ] <- x2
      #new_dm[, wp] <- x2
      
      # Update the energy state: points or point-pairs?
      ppl <- .getPPL(lags = lags, n.lags = n_lags, dist.mat = new_dm, 
                     pairs = pairs)
      new_energy <- .objPPL(n.lags = n_lags, n.pts = n_pts, pairs = pairs,
                            criterion = criterion, distri = distri, ppl = ppl)
                            
      # Evaluate the new system configuration
      random_prob <- ifelse(greedy, 1, stats::runif(1))
      actual_prob <- acceptance$initial * exp(-k / acceptance$cooling)
      if (track) { accept_probs[k] <- actual_prob }
      if (new_energy <= old_energy) {
        old_conf <- new_conf
        old_energy <- new_energy
        count <- 0
        old_dm <- new_dm
      } else {
        if (new_energy > old_energy && random_prob <= actual_prob) {
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
                stopping$max.count, "\n")
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
      if (count == stopping$max.count) {
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
              stopping$max.count, "\n")
        } else {
          break
        }
      }
      if (progress) utils::setTxtProgressBar(pb, k)
     
      if (k == iterations) { break }
    }
    
    # Prepare output
    eval(.prepare_output())
  }
# FUNCTION - CALCULATE THE CRITERION VALUE #####################################
#' @rdname optimPPL
#' @export
objPPL <-
  function (points, candi, 
    # PPL
    lags = 7, lags.type = "exponential", lags.base = 2, cutoff, distri,
    criterion = "distribution", pairs = FALSE,
    # SPSANN
    x.max, x.min, y.max, y.min) {
    
    # Check arguments
    check <- .checkPPL(lags = lags, lags.type = lags.type, pairs = pairs, 
                       lags.base = lags.base, cutoff = cutoff,
                       criterion = criterion, distri = distri, fun = "objPPL")
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare cutoff and lags
    if (missing(cutoff)) {
      eval(.prepare_jittering())
    }
    cutoff <- .cutoffPPL(cutoff = cutoff, x.max = x.max, y.max = y.max)
    lags <- .lagsPPL(lags = lags, lags.type = lags.type, cutoff = cutoff, 
                     lags.base = lags.base)
    n_lags <- length(lags) - 1
    
    # Initial energy state: points or point-pairs
    dm <- SpatialTools::dist1(points[, 2:3])
    ppl <- .getPPL(lags = lags, n.lags = n_lags, dist.mat = dm, pairs = pairs)
    distri <- .distriPPL(n.lags = n_lags, n.pts = n_pts, criterion = criterion,
                         distri = distri, pairs = pairs)
    res <- .objPPL(ppl = ppl, n.lags = n_lags, n.pts = n_pts, pairs = pairs,
                   criterion = criterion, distri = distri)

    # Output
    return (res)
  }
# FUNCTION - POINTS PER LAG-DISTANCE CLASS #####################################
#' @rdname optimPPL
#' @export
countPPL <-
  function (points, candi,
    # PPL
    lags = 7, lags.type = "exponential", lags.base = 2, cutoff, pairs = FALSE,
    # SPSANN
    x.max, x.min, y.max, y.min) {
    
    # Check arguments
    check <- .checkPPL(lags = lags, lags.type = lags.type, pairs = pairs,
                       lags.base = lags.base, cutoff = cutoff, fun = "countPPL")
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare cutoff and lags
    if (missing(cutoff)) {
      eval(.prepare_jittering())
    }
    cutoff <- .cutoffPPL(cutoff = cutoff, x.max = x.max, y.max = y.max)
    lags <- .lagsPPL(lags = lags, lags.type = lags.type, cutoff = cutoff, 
                     lags.base = lags.base)
    n_lags <- length(lags) - 1

    # Distance matrix and counts
    dm <- SpatialTools::dist1(points[, 2:3])
    res <- .getPPL(lags = lags, n.lags = n_lags, dist.mat = dm, pairs = pairs)
    res <- data.frame(lag.lower = lags[-length(lags)], lag.upper = lags[-1],
                      ppl = res)
    return (res)
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
.checkPPL <-
  function (lags, lags.type, lags.base, cutoff, criterion, distri, pairs, fun) {
    
    # Argument 'pairs'
    if (!is.logical(pairs)) {
      res <- c("'pairs' must be a logical value")
      return (res)
    }
    
    # Arguments 'lags' and 'cutoff'
    if (length(lags) == 1) {
#       if (fun != "optimPPL") {
#         if (missing(cutoff)) {
#           res <- c("'cutoff' is mandatory if the lag intervals are not set")
#           return (res)
#         }
#       }
    } else {
      if (!missing(cutoff)) {
        res <- c("'cutoff' cannot be used when the lag intervals are set")
        return (res)
      }
      if (lags[1] == 0) {
        res <- c("lowest lag value must be larger than zero")
        return (res)
      }
    }
    
    # Argument 'lags.type'
    lt <- c("equidistant", "exponential")
    lt <- is.na(any(match(lt, lags.type)))
    if (lt) {
      res <- paste("'lags.type = ", lags.type, "' is not supported", sep = "")
      return (res)
    }
    
    # Argument 'lags.base'
    if (!is.numeric(lags.base) || length(lags.base) > 1) {
      res <- c("'lags.base' must be a numeric value")
      return (res)
    }
    
    if (fun != "countPPL") {
      
      # Argument 'criterion'
      cr <- is.na(any(match(criterion, c("distribution", "minimum"))))
      if (cr) {
        res <- paste("'criterion = ", criterion, "' is not supported", sep = "")
        return (res)
      }
      
      # Argument 'distri'
      if (!missing(distri)) {
        if (!is.numeric(distri)) {
          res <- c("'distri' must be a numeric vector")
          return (res)
        }
        if (length(lags) == 1) {
          if (length(distri) != lags) {
            res <- paste("'distri' must be of length ", lags, sep = "")
            return (res)
          }
        }
        if (length(lags) > 2) {
          nl <- length(lags) - 1
          if (length(distri) != nl) {
            res <- paste("'distri' must be of length ", nl, sep = "")
            return (res)
          }
        }
      }
    }
  }
# INTERNAL FUNCTION - COMPUTE THE DISTRIBUTION #################################
# If required and missing, the wanted distribution of poins (point-pairs) per
# lag distance class is computed internally.
.distriPPL <-
  function (n.lags, n.pts, criterion, distri, pairs) {
    
    if (criterion == "distribution" && missing(distri)) {
     
      if (pairs) { # Point-pairs per lag
        distri <- rep(n.pts * (n.pts - 1) / (2 * n.lags), n.lags)
        
      } else { # Point per lag
        distri <- rep(n.pts, n.lags)
      }
    }
    
    return (distri)
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE ############################
.objPPL <-
  function (ppl, n.lags, n.pts, criterion, distri, pairs) {
    
    if (pairs) { # Point-pairs per lag
      
      if (criterion == "distribution") {
        res <- sum(abs(distri - ppl))
        
      } else { # minimum
        a <- n.pts * (n.pts - 1) / (2 * n.lags)
        res <- a / (min(ppl) + 1)
      }
      
    } else { # Points per lag
      
      if (criterion == "distribution") {
        res <- sum(distri - ppl)
        
      } else { # minimum
        res <- n.pts / (min(ppl) + 1)
      }
    }
    
    # Output
    return (res)
  }
# INTERNAL FUNCTION - NUMBER OF POINTS (POINT-PAIRS) PER LAG-DISTANCE CLASS ####
# Note:
# It is 3 times faster to use the for loop with function 'which()' than using
# 'apply()' with functions 'table()' and 'cut()'.
# apply(X = dist.mat, 1, FUN = function (X) table(cut(X, breaks = lags)))
# apply(X = ppl, 1, FUN = function (X) sum(X != 0))
.getPPL <-
  function (lags, n.lags, dist.mat, pairs) {
    
    ppl <- vector()
    
    if (pairs) { # Point-pairs per lag
      for (i in 1:n.lags) {
        n <- which(dist.mat > lags[i] & dist.mat <= lags[i + 1])
        ppl[i] <- length(n)
      }
      
    } else { # Points per lag
      for (i in 1:n.lags) {
        n <- which(dist.mat > lags[i] & dist.mat <= lags[i + 1], arr.ind = TRUE)
        ppl[i] <- length(unique(c(n)))
      }
    }
    
    return (ppl)
  }
# INTERNAL FUNCTION - BREAKS OF THE LAG-DISTANCE CLASSES #######################
.lagsPPL <-
  function (lags, lags.type, cutoff, lags.base) {
    
    if (length(lags) == 1) {
      
      if (lags.type == "equidistant") {
        lags <- seq(0.0001, cutoff, length.out = lags + 1)
      }
      
      if (lags.type == "exponential") {
        idx <- lags.base ^ c(1:lags - 1)
        lags <- c(0.0001, rev(cutoff / idx))
      }
    }
    
    return (lags)
  }
# INTERNAL FUNCTION - DETERMINE THE CUTOFF #####################################
.cutoffPPL <- 
  function (cutoff, x.max, y.max) {
    
    if (missing(cutoff)) {
      cutoff <- sqrt(x.max ^ 2 + y.max ^ 2)
    }
    
    return (cutoff)
  }
