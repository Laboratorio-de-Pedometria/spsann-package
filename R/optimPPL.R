#' Optimization of sample patterns for for variogram estimation
#' 
#' Optimize a sample pattern for variogram estimation using spatial
#' simulated annealing. The criterion used is the number of
#' points or point-pairs per lag distance class (\code{optimPPL}). 
#' \code{pointsPerLag} and \code{pairsPerLag} count the number of points or 
#' point-pairs per lag distance class. \code{objPoints} and \code{objPairs} 
#' compute the deviation of the observed distribution of counts from a 
#' pre-specified distribution, or the minimum number of points or point pairs 
#' observed over all lag distance classes.
#' 
#' @param lags Integer. The number of lag distance classes. Alternatively, a 
#' vector of numeric values with the lower and upper limits of each lag 
#' distance class. The lowest value must be different from zero, e.g. 0.0001. 
#' Defaults to \code{lags = 7}.
#' 
#' @param lags.type Character value defining the type of lag distance classes. 
#' Available options are \code{"equidistant"}, for equidistant lag distance 
#' classes, and \code{"exponential"}, for exponentially spaced lag distance 
#' classes. Defaults to \code{lags.type = "exponential"}. See \sQuote{Details} 
#' for more information.
#' 
#' @param lags.base Numeric value defining the creation of exponentially spaced
#' lag distance classes. Defaults to \code{lags.base = 2}. See \sQuote{Details}
#' for more information.
#' 
#' @param cutoff Numeric value defining the maximum distance up to which 
#' lag distance classes are created. Used only when lag distance classes are 
#' not defined. See \sQuote{Details} for more information.
#' 
#' @param criterion Character value defining the measure that should be 
#' returned to describe the energy state of the current system configuration. 
#' Available options are \code{"minimum"} and \code{"distribution"}. The first 
#' returns the minimum number of points or point-pairs observed over all lag 
#' distance classes. The second returns the sum of the differences between a 
#' pre-specified distribution and the observed distribution of counts of points
#' or point-pairs per lag distance class. Defaults to 
#' \code{objective = "minimum"}. See \sQuote{Details} for more information.
#' 
#' @param pre.distri Vector of numeric values used to pre-specify the 
#' distribution of points or point-pair with which the observed counts of points
#' or point-pairs per lag distance class is compared. Used only when 
#' \code{criterion = "distribution"}. Defaults to a uniform distribution. See 
#' \sQuote{Details} for more information.
#' 
#' @template spJitter_doc
#' @template spSANN_doc
#' 
#' @details
#' \subsection{Distances}{
#' Euclidean distances between points are calculated. This computation requires
#' the coordinates to be projected. The user is responsible for making sure that
#' this requirement is attained.
#' }
#' \subsection{Distribution}{
#' Using the default uniform distribution means that the desired number of 
#' \strong{point-pairs} per lag distance class is equal to 
#' \eqn{n \times (n - 1) / (2 \times lag)}, where \eqn{n} is the total number
#' of points and \eqn{lag} is the number of lag distance classes.
#' 
#' Using the default uniform distribution means that the desired number of 
#' \strong{points} per lag distance class is equal to the total number of 
#' points. This is the same as expecting that each point contributes 
#' to every lag distance class.
#' 
#' Distributions other than the available options can be easily implemented 
#' changing the arguments \code{lags}, \code{lags.base} and \code{pre.distri}.
#' }
#' \subsection{Type of lags}{
#' Two types of lag distance classes can be created by default. The first 
#' are evenly spaced lags (\code{lags.type = "equidistant"}). They are created 
#' by simply dividing the distance interval from 0.0001 to \code{cutoff} by the
#' required number of lags. The minimum value of 0.0001 guarantees that a point
#' does not form a pair with itself.
#' 
#' The second type of lag distance classes is defined by exponential spacings
#' (\code{lags.type = "exponential"}). The spacings are defined by the base 
#' \eqn{b} of the exponential expression \eqn{b^n}, where \eqn{n} is the 
#' required number of lags. The base is defined using the argument 
#' \code{lags.base}.
#' }
#' \subsection{Criteria}{
#' There are two optimizing criteria implemented. The first is called using
#' \code{criterion = "distribution"} and is used to minimize the sum of
#' differences between a pre-specified distribution and the observed 
#' distribution of points or point-pairs per lag distance class.
#' 
#' The second criterion is called using \code{criterion = "minimum"}. It 
#' corresponds to maximizing the minimum number of points or point-pairs 
#' observed over all lag distance classes.
#' }
#' @return
#' \code{optimPPL} returns a matrix: the optimized sample pattern with 
#' the evolution of the energy state during the optimization as an attribute.
#' 
#' \code{pointsPerLag} and \code{pairsPerLag} return a data.frame with three 
#' columns: a) the lower and b) upper limits of each lag distance class, and 
#' c) the number of points or point-pairs per lag distance class.
#' 
#' \code{objPoints} and \code{objPairs} return a numeric value depending on the
#' choice of \code{criterion}. If \code{criterion = "distribution"}, the sum of
#' the differences between the pre-specified and observed distribution of counts
#' of points or point-pairs per lag distance class. If 
#' \code{criterion = "minimum"}, the inverse of the minimum count of points or 
#' point pairs over all lag distance classes multiplied by a constant.
#' 
#' @references
#' Bresler, E.; Green, R. E. \emph{Soil parameters and sampling scheme for 
#' characterizing soil hydraulic properties of a watershed}. Honolulu: 
#' University of Hawaii at Manoa, p. 42, 1982.
#' 
#' Marler, R. T.; Arora, J. S. Function-transformation methods for 
#' multi-objective optimization. \emph{Engineering Optimization}. v. 37, p. 
#' 551-570, 2005.
#' 
#' Russo, D. Design of an optimal sampling network for estimating the variogram.
#' \emph{Soil Science Society of America Journal}. v. 48, p. 708-716, 1984.
#' 
#' Truong, P. N.; Heuvelink, G. B. M.; Gosling, J. P. Web-based tool for expert
#' elicitation of the variogram. \emph{Computers and Geosciences}. v. 51, p. 
#' 390-399, 2013.
#' 
#' Warrick, A. W.; Myers, D. E. Optimization of sampling locations for variogram
#' calculations. \emph{Water Resources Research}. v. 23, p. 496-500, 1987.
#' 
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @aliases optimPPL pointsPerLag objPoints pairsPerLag objPairs
#' @keywords spatial optimize
#' @concept simulated annealing
#' @export
#' @examples
#' require(sp)
#' data(meuse)
#' meuse <- as.matrix(meuse[, 1:2])
#' meuse <- matrix(cbind(c(1:dim(meuse)[1]), meuse), ncol = 3)
#' pointsPerLag(meuse, cutoff = 1000)
#' objPoints(meuse, cutoff = 1000)
# UNIT TEST ####################################################################
# Use \code{lags = 1} with \code{pointsPerLag} and \code{pairsPerLag} to check
# that the functions are working correctly. They should return the total number 
# of points in \code{points} and the total possible number of point-pairs 
# \eqn{n \times (n - 1) / 2}, respectively.
# FUNCTION - MAIN ##############################################################
optimPPL <-
  function (points, candidates, lags = 7, lags.type = "exponential", 
            lags.base = 2, cutoff = NULL, criterion = "distribution", 
            pre.distri = NULL, x.max, x.min, y.max, y.min, iterations = 10000,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            boundary, progress = TRUE, verbose = TRUE) {
    
    # ASR: check lags... can be a vector with lower and upper limits
    if (ncol(candidates) != 3) stop ("'candidates' must have three columns")
    if (plotit) par0 <- par() # ASR: put on.exit()
    if (is.integer(points)) {
      n_pts <- points
      points <- sample(c(1:dim(candidates)[1]), n_pts)
      points <- candidates[points, ]
    } else {
      n_pts <- nrow(points)
    }
    n_lags <- lags
    lags <- .getLagBreaks(lags, lags.type, cutoff, lags.base)
    sys_config0 <- points
    old_sys_config <- sys_config0
    
    # Initial energy state
    dist_mat <- as.matrix(dist(sys_config0[, 2:3], method = "euclidean"))
    point_per_lag <- .getPointsPerLag(lags, dist_mat)
    energy_state0 <- .objPointsPerLag(point_per_lag, n_lags, n_pts, 
                                      criterion, pre.distri)
    
    # other settings for the simulated annealing algorithm
    old_dist_mat <- dist_mat
    new_dist_mat <- dist_mat
    best_dist_mat <- dist_mat
    count <- 0
    old_energy_state <- energy_state0
    best_energy_state <- Inf
    energy_states <- vector()
    accept_probs <- vector()
    x_max0 <- x.max
    y_max0 <- y.max
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()
    
    # begin the main loop
    for (k in 1:iterations) {
      
      # jitter one of the points and update x.max and y.max
      which_point <- sample(c(1:n_pts), 1)
      new_sys_config <- spJitterFinite(old_sys_config, candidates, x.max, 
                                       x.min, y.max, y.min, which_point)
      x.max <- x_max0 - (k / iterations) * (x_max0 - x.min)
      y.max <- y_max0 - (k / iterations) * (y_max0 - y.min)
      
      # update the distance matrix and calculate the new energy state
      new_dist_mat <- .updatePPLCpp(new_sys_config[, 2:3], old_dist_mat,
                                    which_point)
      point_per_lag <- .getPointsPerLag(lags, new_dist_mat)
      new_energy_state <- .objPointsPerLag(point_per_lag, n_lags, n_pts,
                                           criterion, pre.distri)
      
      # evaluate the new system configuration
      random_prob <- runif(1)
      actual_prob <- acceptance[[1]] * exp(-k / acceptance[[2]])
      accept_probs[k] <- actual_prob
      if (new_energy_state <= old_energy_state) {
        old_sys_config   <- new_sys_config
        old_energy_state <- new_energy_state
        old_dist_mat     <- new_dist_mat
        count <- 0
      } else {
        if (new_energy_state > old_energy_state & random_prob <= actual_prob) {
          old_sys_config   <- new_sys_config
          old_energy_state <- new_energy_state
          old_dist_mat     <- new_dist_mat
          count <- count + 1
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... p = ", 
                random_prob, "\n")
          }
        } else {
          new_energy_state <- old_energy_state
          new_sys_config   <- old_sys_config
          new_dist_mat     <- old_dist_mat
          count <- count + 1
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... stops at",
                stopping[[1]], "\n")
          }
        }
      }
      
      # Best energy state
      energy_states[k] <- new_energy_state
      if (new_energy_state < best_energy_state / 1.0000001) {
        best_k <- k
        best_sys_config       <- new_sys_config
        best_energy_state     <- new_energy_state
        best_old_energy_state <- old_energy_state
        old_sys_config        <- old_sys_config
        best_dist_mat         <- new_dist_mat
        best_old_dist_mat     <- old_dist_mat
      }
      
      # Plotting
      if (plotit && any(round(seq(1, iterations, 10)) == k)) {
        .spSANNplot(energy_state0, energy_states, k, acceptance, 
                    accept_probs, boundary, new_sys_config[, 2:3],
                    sys_config0[, 2:3], y_max0, y.max, x_max0, x.max)
      } 
      
      # Freezing parameters
      if (count == stopping[[1]]) {
        if (new_energy_state > best_energy_state * 1.000001) {
          old_sys_config   <- old_sys_config
          new_sys_config   <- best_sys_config
          new_dist_mat     <- best_dist_mat
          old_energy_state <- best_old_energy_state
          new_energy_state <- best_energy_state
          old_dist_mat     <- best_old_dist_mat
          count <- 0
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
    res <- .spSANNout(new_sys_config, energy_state0, energy_states, time0)
    return (res)
  }
# FUNCTION - CALCULATE THE CRITERION VALUE #####################################
#' @rdname optimPPL
#' @export
objPoints <-
  function (points, lags = 7, lags.type = "exponential", lags.base = 2, 
            cutoff = NULL, criterion = "distribution", pre.distri = NULL) {
    n_pts <- dim(points)[1]
    n_lags <- lags
    lags <- .getLagBreaks(lags, lags.type, cutoff, lags.base)
    dm <- as.matrix(dist(points[, 2:3], method = "euclidean"))
    ppl <- .getPointsPerLag(lags, dm)
    res <- .objPointsPerLag(ppl, n_lags, n_pts, criterion, pre.distri)
    return (res)
  }
# FUNCTION - POINTS PER LAG DISTANCE CLASS #####################################
#' @rdname optimPPL
#' @export
pointsPerLag <-
  function (points, lags = 7, lags.type = "exponential", lags.base = 2, 
            cutoff = NULL) {
    lags <- .getLagBreaks(lags, lags.type, cutoff, lags.base)
    dm <- as.matrix(dist(points[, 2:3], method = "euclidean"))
    res <- .getPointsPerLag(lags, dm)
    res <- data.frame(lag.lower = lags[-length(lags)], points = res, 
                      lag.upper = lags[-1])
    return (res)
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE ############################
.objPointsPerLag <-
  function (points.per.lag, n.lags, n.pts, criterion = "minimum",
            pre.distri = NULL) {
    if (criterion == "distribution") {
      if (!is.null(pre.distri)) {
        if (!is.numeric(pre.distri)) {
          stop ("pre.distri should be of class numeric")
        }
        if (length(pre.distri) != n.lags) {
          stop ("the length of 'pre.distri' should match the number of lags")
        }
      } else {
        pre.distri <- rep(n.pts, n.lags)
      }    
      res <- sum(pre.distri - points.per.lag)
      return (res) 
    }
    if (criterion == "minimum") {
      res <- n.pts / (min(points.per.lag) + 1)
      return (res)
    }
  }
# INTERNAL FUNCTION - NUMBER OF POINTS PER LAG DISTANCE CLASS ##################
# It is 3 times faster to use the for loop with function 'which' than when
# using 'apply' with functions 'table' and 'cut'.
# apply(X = dist.mat, 1, FUN = function (X) table(cut(X, breaks = lags)))
# apply(X = point_per_lag, 1, FUN = function (X) sum(X != 0))
.getPointsPerLag <- function (lags, dist.mat) {
  point_per_lag <- vector()
  for (i in 1:c(length(lags) - 1)) {
    n <- which(dist.mat > lags[i] & dist.mat <= lags[i + 1], arr.ind = TRUE)
    point_per_lag[i] <- length(unique(c(n)))
  }
  return (point_per_lag)
}
# INTERNAL FUNCTION - BREAKS OF THE LAG DISTANCE CLASSES #######################
.getLagBreaks <-
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
# # POINT PAIRS PER LAG DISTANCE CLASS
# .pairsPerLag <- 
#   function (points, lags, lags.type = "equidistant", lags.base = 2,
#             cutoff = NULL) {
#     if (missing(points)) {
#       stop ("'points' is a mandatory argument")
#     }
#     if (missing(lags) || !is.numeric(lags)) {
#       stop ("'lags' should be a numeric value or vector")
#     }
#     if (length(lags) == 1 && is.null(cutoff)) {
#       stop ("'cutoff' is a mandatory when the lag intervals are not specified") 
#     }
#     if (length(lags) > 1 && !is.null(cutoff)) {
#       stop ("'cutoff' cannot be used when the lag intervals are specified")
#     }
#     d <- dist(points, method = "euclidean")
#     if (length(lags) == 1) {
#       if (lags.type == "equidistant") {
#         lags <- seq(0, cutoff, length.out = lags + 1)
#       }
#       if (lags.type == "exponential") {
#         idx <- vector()
#         for (i in 1:lags - 1) {
#           idx[i] <- lags.base ^ i
#         }
#         lags <- c(0, rev(cutoff / idx), cutoff)
#       }
#     }
#     pairs <- vector()
#     for (i in 1:length(lags)) {
#       n <- which(d > lags[i] & d <= lags[i + 1])
#       pairs[i] <- length(n)
#     }
#     res <- data.frame(lag.lower = lags[-length(lags)], 
#                       pairs = pairs[-length(lags)], lag.upper = lags[-1])
#     return (res)
#   }
# # OBJECIVE FUNCTION - POINT PAIRS PER LAG DISTANCE CLASS
# .objPairs <- 
#   function (points, lags, lags.type = "equidistant", lags.base = 2,
#             cutoff = NULL, criterion = "minimum", pre.distri) {
#     if (missing(points)) {
#       stop ("'points' is a mandatory argument")
#     }
#     if (missing(lags) || !is.numeric(lags)) {
#       stop ("'lags' should be a numeric value or vector")
#     }
#     if (length(lags) == 1 && is.null(cutoff)) {
#       stop ("'cutoff' is a mandatory when the lag intervals are not specified") 
#     }
#     if (length(lags) > 1 && !is.null(cutoff)) {
#       stop ("'cutoff' cannot be used when the lag intervals are specified")
#     }
#     n_pts <- dim(points)[1]
#     if (length(lags) > 1) {
#       n_lags <- length(lags) - 1
#     } else {
#       n_lags <- lags
#     }
#     if (criterion == "distribution") {
#       if (!missing(pre.distri)) {
#         if (!is.numeric(pre.distri)) {
#           stop ("pre.distri should be of class numeric")
#         }
#         if (length(pre.distri) != n_lags) {
#           stop ("the length of 'pre.distri' should match the number of lags")
#         }
#       } else {
#         pre.distri <- rep(n_pts * (n_pts - 1) / (2 * n_lags), n_lags)
#       }
#       pairs <- pairsPerLag(points, lags = lags, lags.type = lags.type,
#                            lags.base = lags.base, cutoff = cutoff)
#       res <- sum(pre.distri - pairs$pairs)
#       return (res)
#     }
#     if (criterion == "minimum") {
#       pairs <- pairsPerLag(points, lags = lags, cutoff = cutoff, 
#                            lags.type = lags.type, lags.base = lags.base)
#       res <- 10000 * (min(pairs$pairs) + 1)
#       return (res)
#     }    
#   }
# # End!
