#' Random perturbation of spatial points
#' 
#' Perturb the coordinates of spatial points (\sQuote{jittering}).
#' 
#' @template spJitter_doc
#' 
#' @param points Integer value, integer vector, data frame or matrix. If 
#' \code{points} is an integer value, it defines the number of points that 
#' should be randomly sampled from \code{candi} to form the starting system
#' configuration. If \code{points} is a vector of integer values, it contains 
#' the row indexes of \code{candi} that correspond to the points that form the
#' starting system configuration. If \code{points} is a data frame or matrix, 
#' it must have three columns in the following order: \code{[, "id"]} the 
#' row indexes of \code{candi} that correspond to each point, \code{[, "x"]} 
#' the projected x-coordinates, and \code{[, "y"]} the projected y-coordinates.
#' Note that in the later case, \code{points} must be a subset of \code{candi}.
#'
#' @param candi Data frame or matrix with the candidate locations for the
#' jittered points. \code{candi} must have two columns in the following 
#' order: \code{[, "x"]} the projected x-coordinates, and \code{[, "y"]} the 
#' projected y-coordinates.
#'
#' @param x.max,x.min,y.max,y.min Numeric value defining the minimum and 
#' maximum quantity of random noise to be added to the projected x- and 
#' y-coordinates. The minimum quantity should be equal to, at least, the 
#' minimum distance between two neighbouring candidate locations. The units 
#' are the same as of the projected x- and y-coordinates. If missing, they 
#' are estimated from \code{candi}.
#' 
#' @param cellsize Vector with two numeric values defining the horizontal (x) 
#' and vertical (y) spacing between the candidate locations in \code{candi}. A 
#' single value can be used if the spacing in the x- and y-coordinates is the
#' same.
#' 
#' @param which.point Integer values defining which point should be perturbed. 
#' 
#' @return
#' A matrix with the jittered projected coordinates of the points.
#' 
#' @author 
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' 
#' @seealso \code{ssaOptim}, \code{\link[sp]{zerodist}}, 
#' \code{\link[base]{jitter}}, \code{\link[geoR]{jitter2d}}.
#' @concept jitter perturb
#' @export
#' @examples
#' require(sp)
#' data(meuse.grid)
#' meuse.grid <- as.matrix(meuse.grid[, 1:2])
#' meuse.grid <- matrix(cbind(1:dim(meuse.grid)[1], meuse.grid), ncol = 3)
#' pts1 <- sample(c(1:dim(meuse.grid)[1]), 155)
#' pts2 <- meuse.grid[pts1, ]
#' pts3 <- spJitter(points = pts2, candi = meuse.grid, x.min = 40,
#'                  x.max = 100, y.min = 40, y.max = 100, 
#'                  which.point = 10, cellsize = 40)
#' plot(meuse.grid[, 2:3], asp = 1, pch = 15, col = "gray")
#' points(pts2[, 2:3], col = "red", cex = 0.5)
#' points(pts3[, 2:3], pch = 19, col = "blue", cex = 0.5)
#' 
#' # Cluster of points
#' pts1 <- c(1:55)
#' pts2 <- meuse.grid[pts1, ]
#' pts3 <- spJitter(points = pts2, candi = meuse.grid, x.min = 40,
#'                 x.max = 80, y.min = 40, y.max = 80, 
#'                 which.point = 1, cellsize = 40)
#' plot(meuse.grid[, 2:3], asp = 1, pch = 15, col = "gray")
#' points(pts2[, 2:3], col = "red", cex = 0.5)
#' points(pts3[, 2:3], pch = 19, col = "blue", cex = 0.5)
# FUNCTION #####################################################################
spJitter <-
  function (points, candi, x.max, x.min, y.max, y.min, which.point, cellsize) {
    
    if (length(cellsize) == 1) { cellsize <- rep(cellsize, 2) }
    
    # Get candidate locations using Cpp
    pt1 <- .spJitterCpp(x = points[, 2:3], y = candi[, 2:3], xmax = x.max, 
                        xmin = x.min, ymax = y.max, ymin = y.min, 
                        idx = which.point)
    
    # Get candidate locations
    pt1 <- pt1[pt1 != 0]
    
    # Select one candidate location
    pt2 <- candi[sample(pt1, 1), ]
    
    # Compute the new x- and y-coordinates
    pt3 <- runif(n = 2, min = -0.5, max = 0.5) * cellsize
    pt2[-1] <- pt2[-1] + pt3
    
    # Output
    res <- points
    res[which.point, ] <- pt2
    return (res)
  }
# spJitterFinite ###############################################################
# spJitterFinite <-
#   function (points, candi, x.max, x.min, y.max, y.min, which.point) {
#     
#     # Get candidate locations using Cpp
#     pt1 <- .spJitterCpp(points[, 2:3], candi[, 2:3], x.max, x.min, y.max, 
#                         y.min, which.point)
#       
#     # Get candidate locations
#     pt1 <- pt1[pt1 != 0]
#     
#     # Select one candidate location
#     pt2 <- sample(pt1, 1)
#     
#     # Check if it already is in the sample (duplicated)
#     dup <- duplicated(c(pt2, points[, 1]))
#     
#     # If it already exists, we try to find another point as many times as
#     # there are point in the sample. The reason for this choice is that the 
#     # more points we have, the more likely it is that the candidate point
#     # already is included in the sample.
#     if (any(dup)) {
#       ntry <- 0
#       while (any(dup)) {
#         pt2 <- sample(pt1, 1)
#         dup <- duplicated(c(pt2, points[, 1]))
#         ntry <- ntry + 1
#         if (ntry == length(nrow(points))) {
#           pt2 <- which.point
#           break
#         }
#       }
#     }
#     
#     res <- points
#     res[which.point, ] <- candi[pt2, ]
#     return (res)
#   }
# OLD spJitterFinite FUNCTION ##################################################
# parameters of the old implementation
# \item{where}{An object of class SpatialPolygons defining the spatial domain 
# to which the perturbed points should be constrained. Used only when 
# \code{finite = FALSE}. See \sQuote{Details} for more information.}
# \item{finite}{Logical value defining if the spatial domain is finite or 
# infinite. This is a mandatory argument. See \sQuote{Details} for more 
# information.}
# \item{x.coord}{List with two sub-arguments defining how the x coordinate 
# should be perturbed. \code{min} and \code{max} define the minimum and maximum
# quantity of random noise to be added to the x coordinate. If 
# \code{finite = TRUE}, then \code{min} should be equal to, at least, the 
# minimum distance between two neighboring candidate locations. If 
# \code{finite = FALSE}, then \code{min} should be equal to, at least, the 
# value passed to the argument \code{zero}. This is a mandatory argument. See 
# \sQuote{Details} for more information.}
# \item{y.coord}{The same as for \code{x.coord}.}
# \item{zero}{Numeric value defining the distance less than or equal to which 
# two points are considered to have zero distance. See more information at 
# \code{\link[sp]{zerodist}}. Used only when \code{finite = FALSE}. Defaults 
# to \code{zero = 1}.}
# \item{iterations}{Integer value defining the maximum number of iterations 
# that should be used when constraining the perturbed points to the spatial 
# domain defined by the argument \code{where}. Defaults to 
# \code{iterations = 10000}. Used only when \code{finite = FALSE}. See 
# \sQuote{Details} for more information.}
# \item{verbose}{Logical for printing details about the success of the 
# algorithm when constraining the perturbed points to the spatial domain 
# defined by the argument \code{where}. Used only when 
# \code{finite = FALSE}.}
# \item{size}{Integer value defining the number of points that should be 
# perturbed at each iteration of the simulated annealing exercise 
# (\code{spSANN}). Defaults to \code{size = 1}. See \code{spSANN} for more 
# information.}
# \item{size.factor}{Numeric value defining the factor by which the number of 
# points that should be perturbed at each iteration of the simulated annealing
# exercise (\code{spSANN}) is decreased. Used only when \code{size} is larger
# than 1. See \code{spSANN} for more information.}
# spJitterFinite <-
#   function (points, candi, x.max, x.min, y.max, y.min, which.pts) {
#     d_x <- x.max + x.min
#     d_y <- y.max + y.min
#     pt0 <- points[which.pts, ]
#     d_x <- unlist(c(pt0[1] - d_x, pt0[1] + d_x))
#     d_y <- unlist(c(pt0[2] - d_y, pt0[2] + d_y))
#     pt1 <- which(candi[, 1] >= d_x[1] & candi[, 1] <= d_x[2] &
#                    candi[, 2] >= d_y[1] & candi[, 2] <= d_y[2])
#     pt2 <- candi[sample(pt1, 1), ]
#     dup <- duplicated(rbind(pt2, points))
#     if (any(dup)) {
#       while (any(dup)) {
#         pt2 <- candi[sample(pt1, 1), ]
#         dup <- duplicated(rbind(pt2, points))
#       }
#     }
#     res <- points
#     res[which.pts, ] <- pt2
#     return (res)
#   }
# OLD spJitter FUNCTION ########################################################
# spJitter <- 
#   function (obj, candi = NULL, where = NULL, which = NULL, finite = NULL,
#             x.coord = list(min = NULL, max = NULL), 
#             y.coord = list(min = NULL, max = NULL), 
#             zero = 1, iterations = 10000, verbose = TRUE) {
#     if (is.null(finite)) {
#       stop("'finite' is a mandatory argument")
#     }
#     if (missing(obj)) {
#       stop ("'obj' is a mandatory argument")
#     }
#     if (!is.list(x.coord) || length(x.coord) != 2) {
#       stop ("'x.coord' should be a list with 2 subarguments")
#     } else {
#       if (is.null(x.coord$min) || is.null(x.coord$max)) {
#         stop ("'min' and 'max' are mandatory subarguments for 'x.coord")
#       }
#     }
#     if (!is.list(y.coord) || length(y.coord) != 2) {
#       stop ("'y.coord' should be a list with 2 subarguments")
#     } else {
#       if (is.null(y.coord$min) || is.null(y.coord$max)) {
#         stop ("'min' and 'max' are mandatory subarguments for 'y.coord")
#       }
#     }
#     if (finite) {
#       if (inherits(obj, "SpatialPoints")) {
#         stop ("'obj' should be a vector of indexes")
#       }
#       if (is.null(candi)) {
#         stop ("'candi' is a mandatory argument")
#       }
#       if (!inherits(candi, what = "data.frame")) {
#         stop ("'candi' should be a data.frame")
#       } else {
#         if (dim(candi)[2] != 2 || 
#               any(colnames(candi) != c("x", "y"))) {
#           stop ("'candi' should have 2 columns named 'x' and 'y'")
#         }
#       }
#       if (unique(which == "all")) {
#         stop ("this option is not functional yet")
#         #res <- candi[sample(c(1:length(candi)), length(obj)), ]
#         } else {
#           if (length(which) > 1) {
#             stop ("this option is not functional yet")
#             #pt0 <- coordinates(obj[which, ])
#             #res <- list()
#             #for (i in 1:length(which)) {
#             #  cand <- coordinates(candi)
#             #  d_x <- x.coord$max + x.coord$min
#             #  d_y <- y.coord$max + y.coord$min
#             #  d_x <- c(pt0[i, "x"] - d_x, pt0[i, "x"] + d_x)
#             #  d_y <- c(pt0[i, "y"] - d_y, pt0[i, "y"] + d_y)
#             #  cand <- which(cand[, "x"] >= d_x[1] & cand[, "x"] <= d_x[2] &
#             #                  cand[, "y"] >= d_y[1] & cand[, "y"] <= d_y[2])
#             #  res[i] <- candi[sample(cand, 1), ]
#             #}
#             #res <- data.frame(t(sapply(res, coordinates)))
#             #colnames(res) <- c("x", "y")
#             #coordinates(res) <- ~ x + y
#             #proj4string(res) <- proj4string(obj)
#             #res <- rbind(obj[-which, ], res)
#           } else {
#             res <- .spJitterFiniteOne(obj = obj, candi = candi,
#                                       x.max = x.coord$max, x.min = x.coord$min,
#                                       y.max = y.coord$max, y.min = y.coord$min,
#                                       which = which)
#           }
#         }
#       return (res)
#       } else {
#         if (is.null(where)) {
#           stop ("'where' is a mandatory argument")
#         } else {
#           if (!inherits(where, what = "SpatialPolygons")) {
#             stop ("'where' should be of class SpatialPolygons")
#           } else {
#             where <- as(where, "SpatialPolygons")
#           }
#         }
#         if (!is.numeric(zero)) {
#           stop ("'zero' should be a numeric value")
#         }
#         if (!is.numeric(iterations)) {
#           stop ("'iterations' should be a numeric value")
#         }
#         if (!inherits(obj, "SpatialPoints") || is.na(proj4string(obj)) || 
#               !is.projected(obj)) {
#           stop ("'obj' should be of class SpatialPoints with projected CRS")
#         } else {
#           obj <- as(obj, "SpatialPoints")
#           colnames(obj@coords) <- c("x", "y")
#           rownames(obj@bbox) <- c("x", "y")
#         }
#         if (unique(which == "all")) {
#           x0 <- coordinates(obj)[, "x"]
#           y0 <- coordinates(obj)[, "y"]
#         }
#         if (is.numeric(which)) {
#           x0 <- coordinates(obj)[which, "x"]
#           y0 <- coordinates(obj)[which, "y"]
#         }
#         x1 <- jitter(x = x0, amount = x.coord$max)
#         dx <- abs(x0 - x1)
#         dx <- which(dx < x.coord$min)
#         x1[dx] <- x0[dx] + sign(x0[dx] - x1[dx]) * x.coord$min
#         y1 <- jitter(x = y0, amount = y.coord$max)
#         dy <- abs(y0 - y1)
#         dy <- which(dy < y.coord$min)
#         y1[dy] <- y0[dy] + sign(y0[dy] - y1[dy]) * y.coord$min
#         res <- data.frame(x = x1, y = y1)
#         coordinates(res) <- ~ x + y
#         proj4string(res) <- proj4string(obj)
#         if (is.numeric(which)) {
#           res <- rbind(obj[-which, ], res)
#         }
#         if (!is.null(where)) {
#           out <- which(gContains(where, res, byid = TRUE) == FALSE)
#           n_out <- length(out)
#           n_iter <- 1
#           while (n_out >= 1) {
#             res_out <- obj[out, ]
#             res <- res[-out, ]
#             x0 <- coordinates(res_out)[, "x"]
#             y0 <- coordinates(res_out)[, "y"]
#             x1 <- jitter(x = x0, amount = x.coord$max)
#             dx <- abs(x0 - x1)
#             dx <- which(dx < x.coord$min)
#             x1[dx] <- x0[dx] + sign(x0[dx] - x1[dx]) * x.coord$min
#             y1 <- jitter(x = y0, amount = y.coord$max)
#             dy <- abs(y0 - y1)
#             dy <- which(dy < y.coord$min)
#             y1[dy] <- y0[dy] + sign(y0[dy] - y1[dy]) * y.coord$min
#             res_out <- data.frame(x = x1, y = y1)
#             coordinates(res_out) <- ~ x + y
#             proj4string(res_out) <- proj4string(res)
#             res <- rbind(res, res_out)
#             while(length(zerodist(res, zero = zero))[1] > 0){
#               dup <- unique(c(zerodist(res, zero = zero)))
#               res_dup <- res[dup, ]
#               res <- res[-dup]
#               x0 <- coordinates(res_dup)[, "x"]
#               y0 <- coordinates(res_dup)[, "y"]
#               x1 <- jitter(x = x0, amount = x.coord$min)
#               y1 <- jitter(x = y0, amount = y.coord$min)
#               res_dup <- data.frame(x = x1, y = y1)
#               coordinates(res_dup) <- ~ x + y
#               proj4string(res_dup) <- proj4string(res)
#               res <- rbind(res, res_dup)
#             }
#             out <- which(gContains(where, res, byid = TRUE) == FALSE)
#             n_out <- length(out)
#             n_iter <- n_iter + 1
#             if (n_iter == iterations) {
#               break
#             }
#           }
#         }
#         if (!is.null(where)) {
#           if (n_out >= 1) {
#             message(paste("spJitter DID NOT converge after ", n_iter, 
#                           " iterations", sep = ""))
#             if (n_out > 1) {
#               stop(paste("there are ", n_out, 
#                          " points outside the spatial domain", sep = ""))
#             } else {
#               stop(paste("there is ", n_out, 
#                          " point outside the spatial domain", sep = ""))
#             }          
#           } else {
#             if (verbose) {
#               message(paste("spJitter converged after ", n_iter,
#                             " iterations", sep = ""))
#             }
#           }
#         }
#         return (res)
#       }
#   }
# End!
