#' Optimization of sample patterns for variogram estimation
#'
#' Optimize a sample pattern for variogram estimation. A criterion is defined so
#' that the optimized sample pattern has a given number of points or point-pairs 
#' contributing to each lag distance class.
#' 
#' @template spJitter_doc
#' @template spSANN_doc
#' 
#' @param lags Integer value. The number of lag distance classes. Alternatively,
#' a vector of numeric values with the lower and upper limits of each lag
#' distance class. The lowest value must be larger than zero, e.g. 0.0001.
#' Defaults to \code{lags = 7}.
#' 
#' @param lags.type Character value. The type of lag distance classes, with
#' options \code{"equidistant"} and \code{"exponential"}. Defaults to
#' \code{lags.type = "exponential"}.
#'
#' @param lags.base Numeric value. Base of the exponential expression used to
#' create exponentially spaced lag distance classes. Used only when 
#' \code{lags.type = "exponential"}. Defaults to \code{lags.base = 2}.
#'
#' @param cutoff Numeric value. The maximum distance up to which lag distance
#' classes are created. Used only when \code{lags} is an integer value. 
#'
#' @param criterion Character value. The measure used to describe the
#' energy state of the current system configuration, with options 
#' \code{"minimum"} and \code{"distribution"}. Defaults to 
#' \code{objective = "distribution"}.
#'
#' @param pre.distri Numeric vector. The pre-specified distribution of points
#' or point-pair with which the observed counts of points or point-pairs per
#' lag distance class is compared. Used only when
#' \code{criterion = "distribution"}. Defaults to a uniform distribution.
#'
#' @details
#' \strong{Distances}: Euclidean distances between points are calculated. This
#' computation requires the coordinates to be projected. The user is responsible
#' for making sure that this requirement is attained.
#' 
#' \strong{Distribution}: Using the default uniform distribution means that the
#' number of point-pairs per lag is equal to 
#' \eqn{n \times (n - 1) / (2 \times lag)}, where \eqn{n} is the total number
#' of points and \eqn{lag} is the number of lags. Using the default uniform
#' distribution means that the number of points per lag is equal to the total
#' number of points. This is the same as expecting that each point contributes
#' to every lag. Distributions other than the available options can be easily
#' implemented changing the arguments \code{lags}, \code{lags.base} and
#' \code{pre.distri}.
#' 
#' \strong{Type of lags}: Two types of lag distance classes can be created by
#' default. The first are evenly spaced lags (\code{lags.type = "equidistant"}).
#' They are created by simply dividing the distance interval from 0.0001 to
#' \code{cutoff} by the required number of lags. The minimum value of 0.0001
#' guarantees that a point does not form a pair with itself. The second type of
#' lags is defined by exponential spacings (\code{lags.type = "exponential"}).
#' The spacings are defined by the base \eqn{b} of the exponential expression
#' \eqn{b^n}, where \eqn{n} is the required number of lags. The base is defined
#' using the argument \code{lags.base}.
#'
#' \strong{Criteria}: There are two optimizing criteria implemented. The first
#' is called using \code{criterion = "distribution"} and is used to minimize the
#' sum of differences between a pre-specified distribution and the observed 
#' distribution of points or point-pairs per lag distance class. The second
#' criterion is called using \code{criterion = "minimum"}. It corresponds to
#' maximizing the minimum number of points or point-pairs observed over all lag
#' distance classes.
#' 
#' @return
#' \code{optimPPL} returns a matrix: the optimized sample pattern with
#' the evolution of the energy state during the optimization as an attribute.
#'
#' \code{pointsPerLag} and \code{pairsPerLag} return a data.frame with three
#' columns: a) the lower and b) upper limits of each lag, and c) the number of
#' points or point-pairs per lag.
#'
#' \code{objPoints} and \code{objPairs} return a numeric value depending on the
#' choice of \code{criterion}. If \code{criterion = "distribution"}, the sum of
#' the differences between the pre-specified and observed distribution of counts
#' of points or point-pairs per lag. If \code{criterion = "minimum"}, the
#' inverse of the minimum count of points or point pairs over all lags
#' multiplied by a constant.
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
#' cutoff <- sqrt((x.max * x.max) + (y.max * y.max))
#' points <- 100
#' set.seed(2001)
#' res <- optimPPL(points = points, candi = candi, lags = 7, lags.base = 2,
#'                 criterion = "distribution", lags.type = "exponential",
#'                 cutoff = cutoff, x.max = x.max, x.min = 40, y.max = y.max, 
#'                 y.min = 40, boundary = boundary, iterations = 100, 
#'                 plotit = TRUE)
#' pointsPerLag(points = res, lags = 7, lags.type = "exponential", 
#' lags.base = 2, cutoff = cutoff)
#' tail(attr(res, "energy.state"), 1) # 92
#' objPoints(points = res, lags = 7, lags.type = "exponential", lags.base = 2,
#'           cutoff = cutoff, criterion = "distribution")
# UNIT TEST ####################################################################
# Use \code{lags = 1} with \code{pointsPerLag} and \code{pairsPerLag} to check
# that the functions are working correctly. They should return the total number
# of points in \code{points} and the total possible number of point-pairs
# \eqn{n \times (n - 1) / 2}, respectively.
# FUNCTION - MAIN ##############################################################
optimPPL <-
  function (points, candi, lags = 7, lags.type = "exponential",
            lags.base = 2, cutoff = NULL, criterion = "distribution",
            pre.distri = NULL, x.max, x.min, y.max, y.min, iterations = 10000,
            acceptance = list(initial = 0.99, cooling = iterations / 10),
            stopping = list(max.count = iterations / 10), plotit = TRUE,
            boundary, progress = TRUE, verbose = TRUE) {
    
    # Check arguments
    # http://www.r-bloggers.com/a-warning-about-warning/
    check <- .spSANNcheck(points, candi, x.max, x.min, y.max, y.min,
                          iterations, acceptance, stopping, plotit, boundary,
                          progress, verbose)
    if (!is.null(check)) stop (check, call. = FALSE)
    check <- .optimPPLcheck(lags, lags.type, lags.base, cutoff, criterion,
                            pre.distri)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    if (plotit) {
      par0 <- par()
      on.exit(suppressWarnings(par(par0)))
    }
    
    # Prepare points
    n_candi <- nrow(candi)
    points <- .spsannPoints(points = points, candi = candi, n.candi = n_candi)
    n_pts <- nrow(points)
    
    # Prepare lags
    if (length(lags) >= 3) {
      n_lags <- length(lags) - 1
    } else {
      n_lags <- lags
      lags <- .getLagBreaks(lags = lags, lags.type = lags.type, 
                            cutoff = cutoff, lags.base = lags.base)
    }
    conf0 <- points
    old_conf <- conf0
    
    # Initial energy state
    dm <- SpatialTools::dist1(conf0[, 2:3])
    ppl <- .getPointsPerLag(lags = lags, dist.mat = dm)
    energy0 <- .objPointsPerLag(ppl = ppl, n.lags = n_lags, n.pts = n_pts,
                                criterion = criterion, pre.distri = pre.distri)
    
    # other settings for the simulated annealing algorithm
    old_dm <- dm
    best_dm <- dm
    count <- 0
    old_energy <- energy0
    best_energy <- Inf
    energies <- vector()
    accept_probs <- vector()
    x_max0 <- x.max
    y_max0 <- y.max
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()
    
    # begin the main loop
    for (k in 1:iterations) {
      
      # jitter one of the points and update x.max and y.max
      wp <- sample(1:n_pts, 1)
      new_conf <- spJitterFinite(points = old_conf, candi = candi,
                                 x.max = x.max, x.min = x.min, y.max = y.max,
                                 y.min = y.min, which.point = wp)
      x.max <- x_max0 - (k / iterations) * (x_max0 - x.min)
      y.max <- y_max0 - (k / iterations) * (y_max0 - y.min)
      
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
      
      # Update the energy state
      ppl <- .getPointsPerLag(lags = lags, dist.mat = new_dm)
      new_energy <- .objPointsPerLag(ppl = ppl, n.lags = n_lags, n.pts = n_pts,
                                     criterion = criterion, 
                                     pre.distri = pre.distri)      
      
      # ASR: This is to test the 'update' function
      #a <- objPoints(points = new_conf, lags = lags, criterion = criterion)
      #if (new_energy != a) {
      # print(new_energy)
      # print(a)
      # break
      #}
      
      # Evaluate the new system configuration
      random_prob <- runif(1)
      actual_prob <- acceptance$initial * exp(-k / acceptance$cooling)
      accept_probs[k] <- actual_prob
      if (new_energy <= old_energy) {
        # Always accepts a better energy
        old_conf   <- new_conf
        old_energy <- new_energy
        count      <- 0
        old_dm     <- new_dm
      } else {
        if (new_energy > old_energy && random_prob <= actual_prob) {
          # Accepts a worse energy depending on the probability
          old_conf   <- new_conf
          old_energy <- new_energy
          count      <- count + 1
          old_dm     <- new_dm
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
                stopping$max.count, "\n")
          }
        }
      }
      
      # Best energy state
      energies[k] <- new_energy
      if (new_energy < best_energy / 1.0000001) {
        best_k          <- k
        best_conf       <- new_conf
        best_energy     <- new_energy
        best_old_energy <- old_energy
        old_conf        <- old_conf
        best_dm         <- new_dm
        best_old_dm     <- old_dm
      }

      # Plotting
      #if (plotit && any(round(seq(1, iterations, 10)) == k)) {
      if (plotit && pedometrics::is.numint(k / 10)) {
        .spSANNplot(energy0 = energy0, energies = energies, k = k, 
                    acceptance = acceptance, accept_probs = accept_probs,
                    boundary = boundary, new_conf = new_conf[, 2:3], 
                    conf0 = conf0[, 2:3], y_max0 = y_max0, y.max = y.max,
                    x_max0 = x_max0, x.max = x.max)
      }

      # Freezing parameters
      if (count == stopping$max.count) {
        if (new_energy > best_energy * 1.000001) {
          old_conf    <- old_conf
          new_conf    <- best_conf
          old_energy  <- best_old_energy
          new_energy  <- best_energy
          count       <- 0
          energies[k] <- new_energy
          new_dm      <- best_dm
          old_dm      <- best_old_dm
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
# FUNCTION - CALCULATE THE CRITERION VALUE #####################################
#' @rdname optimPPL
#' @export
objPoints <-
  function (points, candi, lags = 7, lags.type = "exponential",
            lags.base = 2, cutoff = NULL, criterion = "distribution",
            pre.distri = NULL) {

    # Prepare points
    if (is.integer(points) || is.numint(points)) {
      n_candi <- nrow(candi)
      if (length(points) > 1) {
        n_pts <- length(points)
        points <- candi[points, ]
      }
      if (length(points) == 1) {
        n_pts <- points
        points <- sample(1:n_candi, n_pts)
        points <- candi[points, ] 
      }
    } else {
      n_pts <- nrow(points)
    }

    # Prepare lags
    if (length(lags) >= 3) {
      n_lags <- length(lags) - 1
    } else {
      n_lags <- lags
      lags <- .getLagBreaks(lags = lags, lags.type = lags.type, 
                            cutoff = cutoff, lags.base = lags.base)
    }

    dm <- as.matrix(dist(points[, 2:3], method = "euclidean"))
    ppl <- .getPointsPerLag(lags = lags, dist.mat = dm)
    res <- .objPointsPerLag(ppl = ppl, n.lags = n_lags, n.pts = n_pts, 
                            criterion = criterion, pre.distri = pre.distri)
    return (res)
  }
# FUNCTION - POINTS PER LAG DISTANCE CLASS #####################################
#' @rdname optimPPL
#' @export
pointsPerLag <-
  function (points, candi, lags = 7, lags.type = "exponential",
            lags.base = 2, cutoff = NULL) {

    # Prepare points
    if (is.integer(points) || is.numint(points)) {
      n_pts <- points
      n_candi <- nrow(candi)
      points <- sample(1:n_candi, n_pts)
      points <- candi[points, ]
    } else {
      n_pts <- nrow(points)
    }

    # Prepare lags
    if (length(lags) == 1) {
      lags <- .getLagBreaks(lags, lags.type, cutoff, lags.base)
    }

    dm <- as.matrix(dist(points[, 2:3], method = "euclidean"))
    res <- .getPointsPerLag(lags, dm)
    res <- data.frame(lag.lower = lags[-length(lags)], points = res,
                      lag.upper = lags[-1])
    return (res)
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
.optimPPLcheck <-
  function (lags, lags.type, lags.base, cutoff, criterion, pre.distri) {

    # lags and cutoff
    if (length(lags) == 1 && is.null(cutoff)) {
      res <- paste("'cutoff' is mandatory when the lag intervals are not set")
      return (res)
    }
    if (length(lags) > 1 && !is.null(cutoff)) {
      res <- paste("'cutoff' cannot be used when the lag intervals are set")
      return (res)
    }

    # lags.type
    lt <- c("equidistant", "exponential")
    lt <- is.na(any(match(lt, lags.type)))
    if (lt) {
      res <- paste("'lags.type = ", lags.type, "' is not supported", sep = "")
      return (res)
    }

    # lags.base
    if (!is.numeric(lags.base) || length(lags.base) > 1) {
      res <- paste("'lags.base' must be a numeric value")
      return (res)
    }

    # criterion
    cr <- c("distribution", "minimum")
    cr <- is.na(any(match(cr, criterion)))
    if (cr) {
      res <- paste("'criterion = ", criterion, "' is not supported", sep = "")
      return (res)
    }

    # pre.distri
    if (!is.null(pre.distri)) {
      if (!is.numeric(pre.distri)) {
        res <- paste("'pre.distri' must be a numeric vector")
        return (res)
      }
      if (length(lags) == 1) {
        if (length(pre.distri) != lags) {
          res <- paste("'pre.distri' must be of length ", lags, sep = "")
          return (res)
        }
      }
      if (length(lags) > 2) {
        nl <- length(lags) - 1
        if (length(pre.distri) != nl) {
          res <- paste("'pre.distri' must be of length ", nl, sep = "")
          return (res)
        }
      }
    }
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE ############################
.objPointsPerLag <-
  function (ppl, n.lags, n.pts, criterion, pre.distri = NULL) {
    if (criterion == "distribution") {
      if (is.null(pre.distri)) {
        pre.distri <- rep(n.pts, n.lags)
      }
      res <- sum(pre.distri - ppl)
      } else {
        if (criterion == "minimum") {
          res <- n.pts / (min(ppl) + 1)
        }
      }
    return (res)
  }
# INTERNAL FUNCTION - NUMBER OF POINTS PER LAG DISTANCE CLASS ##################
# It is 3 times faster to use the for loop with function 'which' than when
# using 'apply' with functions 'table' and 'cut'.
# apply(X = dist.mat, 1, FUN = function (X) table(cut(X, breaks = lags)))
# apply(X = ppl, 1, FUN = function (X) sum(X != 0))
.getPointsPerLag <- function (lags, dist.mat) {
  ppl <- vector()
  for (i in 1:c(length(lags) - 1)) {
    n <- which(dist.mat > lags[i] & dist.mat <= lags[i + 1], arr.ind = TRUE)
    ppl[i] <- length(unique(c(n)))
  }
  return (ppl)
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
