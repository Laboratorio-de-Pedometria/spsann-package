#' @importFrom sp bbox
# INTERNAL FUNCTION - PREPARE POINTS ###########################################
.spsannPoints <-
  function (points, candi, n.candi) {
    
    if (is.integer(points) || is.numint(points)) {
      
      # Integer vector
      if (length(points) > 1) {
        points <- candi[points, ]
      }
      
      # Integer value
      if (length(points) == 1) {
        points <- sample(1:n.candi, points)
        points <- candi[points, ] 
      }
    } else {
      # Data frame of matrix
      points <- points
    }
    return (points)
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
.spSANNcheck <-
  function (points, candi, x.max, x.min, y.max, y.min, iterations,
            acceptance, stopping, plotit, boundary, progress, verbose) {
    
    # missing arguments
    arg <- c("points", "candi", "x.max", "x.min", "y.max", "y.min")
    mis <- c(missing(points), missing(candi), missing(x.max), 
             missing(x.min), missing(y.max), missing(y.min))
    if (any(mis)) {
      i <- which(mis == TRUE)
      res <- paste("missing argument: ", arg[i], "\n", sep = "")
      return (res)
    }
    
    # candi
    if (ncol(candi) != 3) {
      res <- paste("'candi' must have three columns")
      return (res)
    }
    
    # iterations
    if (!is.numint(iterations) || length(iterations) > 1) {
      res <- paste("'iterations' must be an integer value")
      return (res)
    }
    
    # acceptance
    aa <- !is.list(acceptance)
    bb <- length(acceptance) != 2
    cc <- is.null(names(acceptance))
    dd <- !all(c(names(acceptance) == c("initial", "cooling")) == TRUE)
    if (aa || bb || cc || dd) {
      res <- paste("'acceptance' must be a list with two named sub-arguments: ",
                   "'initial' and 'cooling'", sep = "")
      return (res)
    }
    
    # stopping
    aa <- !is.list(stopping)
    bb <- length(stopping) != 1
    cc <- is.null(names(stopping))
    dd <- !all(c(names(stopping) == "max.count") == TRUE)
    if (aa || bb || cc || dd) {
      res <- paste("'stopping' must be a list with one named sub-argument: ",
                   "'max.count'", sep = "")
      return (res)
    }
    
    # boundary
    if (plotit) {
      if (missing(boundary)) {
        res <- paste("'boundary' is mandatory when 'plotit = TRUE'")
        return (res)
      }
      if (inherits(boundary, "SpatialPolygon")) {
        res <- paste("'boundary' must be a SpatialPolygon")
        return (res)
      }
    }    
  }
# INTERNAL FUNCTION - PLOTTING #################################################
.spSANNplot <-
  function (energy0, energies, k, acceptance, accept_probs, boundary, new_conf,
            conf0, y_max0, y.max, x_max0, x.max, ...) {
    par(mfrow = c(1, 2))
    
    # plot energy states
    a <- c(energy0, energies[1:k])
    plot(a ~ c(0:k), type = "l", xlab = "iteration", ylab = "energy state", ...)
    abline(h = energy0, col = "red")
    
    # plot acceptance probability
    a <- c(acceptance[[1]], accept_probs[1:k])
    par(new = TRUE)
    plot(a ~ c(0:k), type = "l", axes = FALSE, bty = "n", xlab = "", 
         ylab = "", col = "blue", ylim = c(0, acceptance[[1]]))
    axis(side = 4, at = pretty(range(a)))
    mtext("acceptance probability", side = 4, line = 3)
    
    # plot sample points
    bb <- bbox(boundary)
    plot(boundary)
    points(conf0[, 1], conf0[, 2], pch = 1, cex = 0.5, 
           col = "lightgray")
    points(new_conf[, 1], new_conf[, 2], pch = 20, cex = 0.5)
    
    # plot maximum shift in the x and y coordinates
    x <- c(bb[1, 1], bb[1, 2])
    y <- rep(bb[2, 1], 2) - 0.02 * y_max0
    lines(x = x, y = y, col = "gray", lwd = 12)
    
    y <- c(bb[2, 1], bb[2, 2])
    x <- rep(bb[1, 1], 2) - 0.02 * x_max0
    lines(y = y, x = x, col = "gray", lwd = 12)
    
    x <- c(bb[1, 1], bb[1, 1] + x.max)
    y <- rep(bb[2, 1], 2) - 0.02 * y_max0
    lines(x = x, y = y, col = "orange", lwd = 12)
    
    x <- rep(bb[1, 1], 2) - 0.02 * x_max0
    y <- c(bb[2, 1], bb[2, 1] + y.max)
    lines(y = y, x = x, col = "orange", lwd = 12)
    
    # plot labels for maximum shift in the x and y coordinates
    x <- bb[1, 1] + (bb[1, 2] - bb[1, 1]) / 2
    y <- bb[2, 1] - 0.02 * y_max0
    text(x = x, y = y, labels = "maximum shift in the X axis")
    
    x <- bb[1, 1] - 0.02 * x_max0
    y <- bb[2, 1] + (bb[2, 2] - bb[2, 1]) / 2
    text(y = y, x = x, srt = 90, labels = "maximum shift in the Y axis")
  }
# INTERNAL FUNCTION - PREPARE RESULTS ##########################################
.spSANNout <-
  function (new_conf, energy0, energies, time0, nadir) {
    res <- new_conf
    criterion <- c(energy0, energies)
    
    # Prepare attributes: energy states and running time
    a <- attributes(res)
    a$energy.state <- criterion
    a$iterations <- length(energies)
    running_time <- (proc.time() - time0) / 60
    a$running.time <- running_time
    if (!missing(nadir)) {
      a$nadir <- nadir
    }
    attributes(res) <- a
    
    # Print the number of iterations and running time
    cat("iterations = ", a$iterations, "\n", sep = "")
    cat("running time = ", round(running_time[3], 2), " minutes", sep = "")
    return (res)
  }
# THE ORIGINAL spSANN FUNCTION #################################################
# .energyState <- 
#   function (fun, points, ...) {
#     if (missing(fun) || missing(points)) {
#       stop ("'fun' and 'points' are mandatory arguments")
#     }
#     return (do.call(fun, list(points, ...)))
#   }
# spatial simulated annealing
# spSANN <-
#   function (points, candi, x.max, x.min, y.max, y.min, fun, ...,
#             iterations = 10000, plotit = TRUE, boundary,
#             acceptance = list(initial = 0.99, cooling = iterations / 10),
#             stopping = list(max.count = 200), progress = TRUE, 
#             verbose = TRUE) {
#     if (plotit){
#       par0 <- par()
#     }
#     n_pts             <- dim(points)[1]
#     conf0       <- points
#     old_sys_config    <- conf0
#     energy0     <- .energyState(fun = fun, points = old_sys_config, ...)
#     old_energy_state  <- energy0
#     count             <- 0
#     best_energy_state <- Inf
#     energies     <- vector()
#     accept_probs      <- vector()
#     x_max             <- vector()
#     y_max             <- vector()
#     x_max0            <- x.max
#     y_max0            <- y.max
#     if (progress) {
#       pb <- txtProgressBar(min = 1, max = iterations, style = 3)
#     }
#     time0             <- proc.time()
#     for (k in 1:iterations) {
#       id <- sample(c(1:n_pts), 1)
#       new_conf <- spJitterFinite(old_sys_config, candi = candi,
#                                        x.max = x.max, x.min = x.min, 
#                                        y.max = y.max, y.min = y.min,
#                                        which.pts = id)
#       x.max <- x_max0 - (k / iterations) * (x_max0 - x.min)
#       x_max[k] <- x.max
#       y.max <- y_max0 - (k / iterations) * (y_max0 - y.min)
#       y_max[k] <- y.max
#       new_energy_state <- .energyState(fun=fun, points = new_conf, ...)
#       random_prob <- runif(1)
#       actual_prob <- acceptance[[1]] * exp(-k / acceptance[[2]])
#       accept_probs[k] <- actual_prob
#       if (new_energy_state <= old_energy_state) {
#         old_sys_config <- new_conf
#         old_energy_state <- new_energy_state
#         count <- 0
#       } else {
#         if (new_energy_state > old_energy_state & random_prob<=actual_prob) {
#           old_sys_config <- new_conf
#           old_energy_state <- new_energy_state
#           count <- count + 1
#           if (verbose) {
#             if (count == 1) {
#               cat("\n", count, "iteration with no improvement... p = ", 
#                   random_prob, "\n")
#             } else {
#               cat("\n", count, "iterations with no improvement... p = ", 
#                   random_prob, "\n")
#             }
#           }
#         } else {
#           new_energy_state <- old_energy_state
#           new_conf <- old_sys_config
#           count <- count + 1
#           if (verbose) {
#             if (count == 1) {
#               cat("\n", count, "iteration with no improvement... stops at",
#                   stopping$max.count, "\n")
#             } else {
#               cat("\n", count, "iterations with no improvement... stops at",
#                   stopping$max.count, "\n")
#             }
#           }
#         }
#       }
#       energies[k] <- new_energy_state
#       if (new_energy_state < best_energy_state / 1.0000001) {
#         best_k <- k
#         best_sys_config <- new_conf
#         best_energy_state <- new_energy_state
#         best_old_energy_state <- old_energy_state
#         old_sys_config <- old_sys_config
#       }
#       if (any(round(seq(1, iterations, 10)) == k)) {
#         if (plotit){
#           .spSANNplot(energy0, energies,k,acceptance,accept_probs, 
#                       boundary, new_conf, conf0, y_max0, y_max, 
#                       x_max0, x_max)
#         } 
#       }
#       if (count == stopping[[1]]) {
#         if (new_energy_state > best_energy_state * 1.000001) {
#           old_sys_config <- old_sys_config
#           new_conf <- best_sys_config
#           old_energy_state <- best_old_energy_state
#           new_energy_state <- best_energy_state
#           count <- 0
#           cat("\n", "reached maximum count with suboptimal system configuration\n")
#           cat("\n", "restarting with previously best system configuration\n")
#           if (count == 1) {
#             cat("\n", count, "iteration with no improvement... stops at",
#                 stopping[[1]], "\n")
#           } else {
#             cat("\n", count, "iterations with no improvement... stops at",
#                 stopping[[1]], "\n")
#           }
#         } else {
#           break
#         }
#       }
#       if (progress) {
#         setTxtProgressBar(pb, k)
#       }      
#     }
#     if (progress) {
#       close(pb)
#     }
#     if (plotit){
#       par(par0)
#     }
#     res <- new_conf
#     criterion <- c(energy0, energies)
#     a <- attributes(res)
#     a$energy.state <- criterion
#     running_time <- (proc.time() - time0) / 60
#     a$running.time <- running_time
#     attributes(res) <- a
#     cat("running time = ", round(running_time[3], 2), " minutes", sep = "")
#     return (res)
#   }
# End!
