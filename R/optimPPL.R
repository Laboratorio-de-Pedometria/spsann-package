#' Optimization of sample configurations for variogram identification and estimation
#'
#' Optimize a sample configuration for variogram identification and estimation. A criterion is defined so that
#' the optimized sample configuration has a given number of points or point-pairs contributing to each 
#' lag-distance class (\bold{PPL}).
#' 
#' @inheritParams spJitter
#' @template spSANN_doc
#' @template PPL_doc
#' @template spJitter_doc
#' 
#' @details 
#' \subsection{Lag-distance classes}{
#' Two types of lag-distance classes can be created by default. The first are evenly spaced lags 
#' (\code{lags.type = "equidistant"}). They are created by simply dividing the distance interval from 0.0001 
#' to \code{cutoff} by the required number of lags. The minimum value of 0.0001 guarantees that a point does 
#' not form a pair with itself. The second type of lags is defined by exponential spacings 
#' (\code{lags.type = "exponential"}). The spacings are defined by the base \eqn{b} of the exponential 
#' expression \eqn{b^n}, where \eqn{n} is the required number of lags. The base is defined using the argument 
#' \code{lags.base}. See \code{\link[pedometrics]{vgmLags}} for other details.
#'
#' Using the default uniform distribution means that the number of point-pairs per lag-distance class 
#' (\code{pairs = TRUE}) is equal to \eqn{n \times (n - 1) / (2 \times lag)}, where \eqn{n} is the total 
#' number of points and \eqn{lag} is the number of lags. If \code{pairs = FALSE}, then it means that the 
#' number of points per lag is equal to the total number of points. This is the same as expecting that each 
#' point contributes to every lag. Distributions other than the available options can be easily implemented 
#' changing the arguments \code{lags} and \code{distri}.
#' 
#' There are two optimizing criteria implemented. The first is called using \code{criterion = "distribution"}
#' and is used to minimize the sum of the absolute differences between a pre-specified distribution and the
#' observed distribution of points or point-pairs per lag-distance class. The second criterion is called using
#' \code{criterion = "minimum"}. It corresponds to maximizing the minimum number of points or point-pairs 
#' observed over all lag-distance classes.
#' }
#' 
#' @return
#' \code{optimPPL} returns an object of class \code{OptimizedSampleConfiguration}: the optimized sample
#' configuration with details about the optimization.
#'
#' \code{objPPL} returns a numeric value: the energy state of the sample configuration -- the objective 
#' function value.
#'
#' \code{countPPL} returns a data.frame with three columns: a) the lower and b) upper limits of each 
#' lag-distance class, and c) the number of points or point-pairs per lag-distance class.
#'
#' @references
#' Bresler, E.; Green, R. E. \emph{Soil parameters and sampling scheme for characterizing soil hydraulic
#' properties of a watershed}. Honolulu: University of Hawaii at Manoa, p. 42, 1982.
#'
#' Pettitt, A. N.; McBratney, A. B. Sampling designs for estimating spatial variance components. 
#' \emph{Applied Statistics}. v. 42, p. 185, 1993.
#'
#' Russo, D. Design of an optimal sampling network for estimating the variogram. \emph{Soil Science Society of
#' America Journal}. v. 48, p. 708-716, 1984.
#'
#' Truong, P. N.; Heuvelink, G. B. M.; Gosling, J. P. Web-based tool for expert elicitation of the variogram.
#' \emph{Computers and Geosciences}. v. 51, p.
#' 390-399, 2013.
#'
#' Warrick, A. W.; Myers, D. E. Optimization of sampling locations for variogram calculations. \emph{Water 
#' Resources Research}. v. 23, p. 496-500, 1987.
#'
#' @author
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @aliases optimPPL countPPL objPPL PPL
#' @export
#' @examples
#' #####################################################################
#' # NOTE: The settings below are unlikely to meet your needs.         #
#' #####################################################################
#' \dontrun{
#' # This example takes more than 5 seconds
#' data(meuse.grid, package = "sp")
#' candi <- meuse.grid[, 1:2]
#' schedule <- scheduleSPSANN(
#'   chains = 1,
#'   initial.acceptance = c(0.8, 0.99),
#'   initial.temperature = 9.5,
#'   x.max = 1540, y.max = 2060, x.min = 0,
#'   y.min = 0, cellsize = 40)
#' set.seed(2001)
#' res <- optimPPL(points = 10, candi = candi, schedule = schedule)
#' objSPSANN(res)
#' }
# FUNCTION - MAIN #############################################################################################
optimPPL <-
  function(
    points, candi, lags = 7, lags.type = "exponential", lags.base = 2, cutoff, distri,
    criterion = "distribution", pairs = FALSE, schedule,
    plotit = FALSE, track = FALSE, boundary, progress = "txt", verbose = FALSE) {

    # Objective function name
    objective <- "PPL"

    # Check spsann arguments
    eval(.check_spsann_arguments())

    # Check other arguments
    check <- do.call(.check_ppl_arguments, as.list(match.call()[-1]))
    if (!is.null(check)) {
      stop(check, call. = FALSE)
    }
    # Set plotting options
    eval(.plotting_options())

    # Prepare points and candi
    eval(.prepare_points())

    # Prepare for jittering
    eval(.prepare_jittering())
    
    # compute lags (default behavior)
    lags <- do.call(.compute_variogram_lags, as.list(match.call()[-1]))
    n_lags <- length(lags) - 1

    # Initial energy state: points or point-pairs
    # Use 'old_conf' instead of 'conf0' because the former has information on any existing fixed points.
    # Use 'n_pts + n_fixed_pts' to account for existing fixed points.
    # 
    # dm <- SpatialTools::dist1(conf0[, 2:3])
    # ppl <- .getPPL(lags = lags, n.lags = n_lags, dist.mat = dm, pairs = pairs, n.pts = n_pts)
    # distri <- .distriPPL(
    # n.lags = n_lags, n.pts = n_pts, criterion = criterion, distri = distri, pairs = pairs)
    # energy0 <- data.frame(
    # obj = .objPPL(
    # ppl = ppl, n.lags = n_lags, n.pts = n_pts, criterion = criterion, distri = distri, pairs = pairs))
    dm <- SpatialTools::dist1(old_conf[, 2:3])
    ppl <- .getPPL(lags = lags, n.lags = n_lags, dist.mat = dm, pairs = pairs, n.pts = n_pts + n_fixed_pts)
    distri <- .distriPPL(
      n.lags = n_lags, n.pts = n_pts + n_fixed_pts, criterion = criterion, distri = distri, pairs = pairs)
    energy0 <- data.frame(
      obj = .objPPL(ppl = ppl, n.lags = n_lags, n.pts = n_pts + n_fixed_pts, criterion = criterion, 
                    distri = distri, pairs = pairs))
    
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
          # Use 'n_pts + n_fixed_pts' to account for existing fixed points.
          # 
          # ppl <- .getPPL(
          # lags = lags, n.lags = n_lags, dist.mat = new_dm, pairs = pairs, n.pts = n_pts)
          # new_energy <- data.frame(
          # obj = .objPPL(
          # n.lags = n_lags, n.pts = n_pts, pairs = pairs, criterion = criterion, distri = distri, 
          # ppl = ppl))
          ppl <- .getPPL(
            lags = lags, n.lags = n_lags, dist.mat = new_dm, pairs = pairs, n.pts = n_pts + n_fixed_pts)
          new_energy <- data.frame(
            obj = .objPPL(
              n.lags = n_lags, n.pts = n_pts + n_fixed_pts, pairs = pairs, criterion = criterion, 
              distri = distri, ppl = ppl))
          
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
# FUNCTION - CALCULATE THE OBJECTIVE FUNCTION VALUE ###########################################################
#' @rdname optimPPL
#' @export
objPPL <-
  function(
    points, candi, lags = 7, lags.type = "exponential", lags.base = 2, cutoff, distri, 
    criterion = "distribution", pairs = FALSE, x.max, x.min, y.max, y.min) {
    
    # Check arguments
    check <- do.call(.check_ppl_arguments, as.list(match.call()[-1]))
    if (!is.null(check)) {
      stop(check, call. = FALSE)
    }
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # compute lags (default behavior)
    lags <- do.call(.compute_variogram_lags, as.list(match.call()[-1]))
    n_lags <- length(lags) - 1
    
    # Initial energy state: points or point-pairs
    dm <- SpatialTools::dist1(points[, 2:3])
    ppl <- .getPPL(lags = lags, n.lags = n_lags, dist.mat = dm, pairs = pairs, n.pts = n_pts)
    distri <- .distriPPL(
      n.lags = n_lags, n.pts = n_pts, criterion = criterion, distri = distri, pairs = pairs)
    res <- .objPPL(
      ppl = ppl, n.lags = n_lags, n.pts = n_pts, pairs = pairs, criterion = criterion, distri = distri)
    
    # Output
    return(res)
  }
# FUNCTION - POINTS PER LAG-DISTANCE CLASS ####################################################################
#' @rdname optimPPL
#' @export
countPPL <-
  function(
    points, candi, lags = 7, lags.type = "exponential", lags.base = 2, cutoff, pairs = FALSE, x.max, x.min, 
    y.max, y.min) {

    # Check arguments
    check <- do.call(.check_ppl_arguments, as.list(match.call()[-1]))
    if (!is.null(check)) {
      stop(check, call. = FALSE)
    }
    # Prepare points and candi
    eval(.prepare_points())

    # compute lags (default behavior)
    lags <- do.call(.compute_variogram_lags, as.list(match.call()[-1]))
    n_lags <- length(lags) - 1
    
    # Distance matrix and counts
    dm <- SpatialTools::dist1(points[, 2:3])
    res <- .getPPL(lags = lags, n.lags = n_lags, dist.mat = dm, pairs = pairs, n.pts = n_pts)
    res <- data.frame(lag.lower = lags[-length(lags)], lag.upper = lags[-1], ppl = res)
    
    # Output
    return(res)
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS #########################################################################
# Important: (1) inform default values and (2) allow additional unused arguments
.check_ppl_arguments <-
  function(
    lags = 7, lags.type = "exponential", lags.base = 2, cutoff, criterion = "distribution", distri,
    pairs = FALSE, ...) {
    
    # pairs
    if (!is.logical(pairs)) {
      return(paste0("pairs = '", pairs, "' is not supported"))
    }
    
    # lags and cutoff
    if (length(lags) > 1) {
      if (!missing(cutoff)) {
        return("the cutoff cannot be set when the lag intervals are set too")
      }
      if (lags[1] == 0) {
        return("the lowest lag value must be larger than zero")
      }
    }
    
    # lags.type
    if (!lags.type %in% c("equidistant", "exponential")) {
      return(paste0("lags.type = '", lags.type, "' is not supported"))
    }
    
    # lags.base
    if (!is.numeric(lags.base) || length(lags.base) > 1) {
      return("lags.base must be a numeric value")
    }
    
    # criterion
    if (!criterion %in% c("distribution", "minimum")) {
      return(paste0("criterion = '", criterion, "' is not supported"))
    }
    
    # distri
    if (!missing(distri)) {
      if (!is.numeric(distri)) {
        return("distri must be a vector of numeric values")
      }
      if (length(lags) == 1 & length(distri) != lags) {
        return(paste0("distri must be a vector of length ", lags))
      }
      if (length(lags) > 2 & length(distri) != (length(lags) - 1)) {
        return(paste0("distri must be a vector of length ", (length(lags) - 1)))
      }
    }
  }
# INTERNAL FUNCTION - COMPUTE THE DISTRIBUTION ################################################################
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
  function (lags, n.lags, dist.mat, pairs, n.pts) {
    
    # ppl <- vector()
    
    if (pairs) { # Point-pairs per lag
      # for (i in 1:n.lags) {
      # n <- which(dist.mat > lags[i] & dist.mat <= lags[i + 1])
      # ppl[i] <- length(n)
      # }
      ppl <- diff(sapply(1:length(lags), function (i) length(which(dist.mat <= lags[i]))) - n.pts) * 0.5
    } else { # Points per lag
      # for (i in 1:n.lags) {
      # n <- which(
      # dist.mat > lags[i] & dist.mat <= lags[i + 1], arr.ind = TRUE)
      # ppl[i] <- length(unique(c(n)))
      # }
      ppl <- sapply(1:n.lags, function (i)
        length(unique(c(which(dist.mat > lags[i] & dist.mat <= lags[i + 1], arr.ind = TRUE)))))
    }
    
    return (ppl)
  }
# INTERNAL FUNCTION - BREAKS OF THE LAG-DISTANCE CLASSES ###########################################
.compute_variogram_lags <-
  function(lags = 7, lags.type = "exponential", cutoff, lags.base = 2, x.max, y.max, candi, ...) {
    if (length(lags) == 1) {
      if (missing(cutoff)) {
        if (missing(x.max) & missing(y.max)) {
          x.max <- diff(range(candi[, "x"])) / 2
          y.max <- diff(range(candi[, "y"])) / 2
        }
        cutoff <- sqrt(x.max ^ 2 + y.max ^ 2)
      }
      switch(
        lags.type,
        equidistant = {
          lags <- seq(0.0001, cutoff, length.out = lags + 1)
        },
        exponential = {
          idx <- lags.base ^ c(1:lags - 1)
          lags <- c(0.0001, rev(cutoff / idx))
        })
    }
    return(lags)
  }
