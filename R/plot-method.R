#' Plot an optimized spatial sample configuration
#'
#' Plot the evolution of the objective function value and the optimized spatial
#' sample configuration
#'
# @inheritParams optimACDC
#'
#' @param x Object of class \code{OptimizedSampleConfiguration} returned by one
#' of the \code{optim}-functions.
#'
#' @param which Which plot should be produced: evolution of the objective
#' function value (1), optimized spatial sample configuration (2), or both
#' (1:2)? Defaults to \code{which = 1:2}.
#'
#' @param boundary (Optional) An object of class SpatialPolygons (see sp::SpatialPolygons()) with
#' the outer and inner limits of the spatial sampling domain (see `candi`).
#'
#' @param ... Other options passed to \code{plot}.
#'
#' @rdname plot-method
#' @export
#' @method plot OptimizedSampleConfiguration
#' @aliases plot plot.OptimizedSampleConfiguration
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
#' plot(res)
# MAIN FUNCTION - PLOT OSC #####################################################
plot.OptimizedSampleConfiguration <-
  function(x, which = 1:2, boundary, ...) {

    # Do not try to plot the energy states if they have not been tracked
    if (nrow(x$objective$energy) == 2) {
      which <- 2
    }

    par0 <- graphics::par()
    on.exit(suppressWarnings(graphics::par(par0)))
    if (all(which == 1:2)) {
      graphics::par(mfrow = c(1, 2))
    }

    # Plot the energy states
    if (all(which == 1:2)) {
      # k <- x$spsann$chains[2:3]
      # k <- as.numeric(k[1] * k[2] * nrow(x$points))
      a <- x$objective$energy
      # if (nrow(a) < k) {
        k <- nrow(a) - 1
      # }

      l <- colnames(a)
      n <- ncol(a)
      # Palete based on ColorBrewer (qualitative, Paired, colorblind safe)
      # https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=4
      col <- c(
        "#000000", "#1F78B4", "#33A02C", "#E31A1C", "#A6CEE3",
        "#B2DF8A", "#FB9A99")[1:n]
      if (n > 2) {
        ylim <- range(sapply(a, max))
      } else {
        ylim <- range(a)
      }
      # ylim <- range(sapply(a, max))
      graphics::plot(
        1, type = "n", xlim = c(0, k),
        # ylim = c(0, max(sapply(a, max)) * 1.1),
        ylim = ylim, xlab = "Spatial jitter", ylab = "Objective Function Value")
      graphics::legend(
        "topright", legend = l, lwd = 1, lty = rep(1, n), col = col)
      for (i in 1:ncol(a)) {
        graphics::lines(a[, i] ~ c(0:k), type = "l", col = col[i])
      }
    }

    # Plot optimized sample configuration
    if (all(which == 1:2) || which == 2) {
    # if (which == 1:2 || which == 2) {
      if (!missing(boundary)) {
        bb <- sp::bbox(boundary)
        if (all(x$points$free == 1)) {
          xlab <- ""
        } else {
          xlab <- paste(
            intToUtf8(9679), " = free; ", intToUtf8(10005), " = fixed")
        }
        if (methods::is(boundary, "SpatialPoints")) {
          sp::plot(x = boundary, pch = 20, cex = 0.1, xlab = xlab)
        } else {
          sp::plot(x = boundary, xlab = xlab)
        }
        graphics::points(
          x$points[, "x"], x$points[, "y"],
          pch = ifelse(x$points[, "free"] == 1, 20, 4), cex = 0.75)
      } else {
        graphics::plot(
          x$points[, c("x", "y")],
          pch = ifelse(x$points[, "free"] == 1, 20, 4), cex = 0.75, asp = 1,
          xlab = c("x = fixed, o = free"))
      }
    }
  }
