#' Plot an optimized sample configuration
#' 
#' Plot the evolution of the energy state and the optimized sample configuration
#' 
#' @inheritParams optimACDC
#' 
#' @param osc Object of class \code{'OptimizedSampleConfiguration'} returned by
#' one of the \code{optim}-functions.
#' 
#' @param which Which plot should be produced: evolution of the energy state 
#' (1), optimized sample configuration (2), or both (1:2)? Defaults to
#' \code{which = 1:2}.
#' 
#' @param boundary Object of class Spatial defining the boundary of the 
#' sampling region.
#' 
#' @rdname plot-method
#' @export
#' @examples 
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' covars <- meuse.grid[, 5]
#' schedule <- scheduleSPSANN(initial.temperature = 5, chains = 1,
#'                            x.max = 1540, y.max = 2060, x.min = 0, 
#'                            y.min = 0, cellsize = 40)
#' set.seed(2001)
#' res <- optimCORR(points = 10, candi = candi, covars = covars, 
#'                  use.coords = TRUE, schedule = schedule)
#' plotOSC(res)
# MAIN FUNCTION - PLOT OSC #####################################################
plotOSC <-
  function (osc, which = 1:2, boundary) {
    
    # Do not try to plot the energy states if they have not been traked
    if (nrow(methods::slot(osc, "objective")$energy) == 2) which <- 2
    
    par0 <- graphics::par()
    on.exit(suppressWarnings(graphics::par(par0)))
    if (all(which == 1:2)) {
      graphics::par(mfrow = c(1, 2))
    }
    
    # Plot the energy states
    if (all(which == 1:2)) {
      k <- methods::slot(osc, "spsann")$chains[2:3]
      k <- as.numeric(k[1] * k[2] * nrow(methods::slot(osc, "points")))
      a <- methods::slot(osc, "objective")$energy
      l <- colnames(a)
      n <- ncol(a)
      col <- c("red", rep("black", n - 1))
      ylim <- range(sapply(a, max))
      graphics::plot(1, type = 'n', xlim = c(0, k), 
                     # ylim = c(0, max(sapply(a, max)) * 1.1), 
                     ylim = ylim, 
                     xlab = "jitter", ylab = "energy state")
      graphics::legend("topright", legend = l, lwd = 1, lty = 1:n,
                       col = col)
      for(i in 1:ncol(a)) {
        graphics::lines(a[, i] ~ c(0:k), type = "l", lty = i, 
                        col = col[i])
      }
    }
    
    # Plot optimized sample configuration
    if (which == 1:2 || which == 2) {
      if (!missing(boundary)) {
        bb <- sp::bbox(boundary)
        if (methods::is(boundary, "SpatialPoints")) {
          sp::plot(x = boundary, pch = 20, cex = 0.1)
        } else {
          sp::plot(x = boundary)
        }
        graphics::points(methods::slot(osc, "points")[, "x"], 
                         methods::slot(osc, "points")[, "y"],
                         pch = 20, cex = 0.5)  
      } else {
        graphics::plot(methods::slot(osc, "points")[, c("x", "y")], pch = 20,
                       cex = 0.5, asp = 1)
      }
    }
  }
