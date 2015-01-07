#' @title Random perturbation of spatial points
#' 
#' @description
#' This function perturbs the coordinates of spatial points. This is also known
#' as \sQuote{jittering}.
#' 
#' @template spJitter_doc
#' @param which.point Integer values defining which point should be perturbed. 
#' The current version accepts only one point to be perturbed at a time. See 
#' \sQuote{Details} for more information.
#'
#' @details
#' This function perturbs the coordinates of spatial points adding random noise,
#' a process also known as \sQuote{jittering}. There are two ways of jittering 
#' the coordinates. They differ on how the the set of candidate locations is 
#' defined.
#' 
#' \subsection{Finite set of candidate locations}{
#' The first method uses a finite set of candidate locations for the perturbed 
#' points. This method usually is the fastest because it does not require the 
#' use of complex routines to check if the perturbed point falls inside the 
#' spatial domain. Since the candidate locations is a finite set, any perturbed 
#' point will inexorably fall inside the spatial domain. This is a very 
#' important feature in optimization exercises with complex objective functions
#' such as simulated annealing when repetitive perturbation is required.
#' 
#' The arguments \code{x.min}, \code{y.min}, \code{x.max}, and \code{y.max} are
#' used to define a rectangular window containing the set of effective candidate
#' locations for the point defined with the argument \code{which.point}. The new
#' location is then randomly sampled from the set of effective candidate 
#' locations and checked against existing points to avoid duplicates. The
#' current implementation does not enable to define the direction of the 
#' perturbation, nor to perturb more than one point at a time.
#' }
#' \subsection{Infinite set of candidate locations}{
#' The current version does not accept using an infinite set of candidate 
#' locations.
#' }
#' @return A matrix with the jittered coordinates of the points.
#' 
#' @references
#' Edzer Pebesma, Jon Skoien with contributions from Olivier Baume, A. Chorti, 
#' D.T. Hristopulos, S.J. Melles and G. Spiliopoulos (2013). 
#' \emph{intamapInteractive: procedures for automated interpolation - methods 
#' only to be used interactively, not included in intamap package.} R package 
#' version 1.1-10. \url{http://CRAN.R-project.org/package=intamapInteractive}
#' 
#' Van Groenigen, J.-W. \emph{Constrained optimization of spatial sampling: 
#' a geostatistical approach.} Wageningen: Wageningen University, p. 148, 1999.
#' 
#' @author 
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' 
#' Gerard Heuvelink \email{gerard.heuvelink@@wur.nl}
#' 
#' @note
#' This function is under construction! The current version does not accept 
#' using an infinite set of candidate locations or to perturb 
#' more than one point at a time.
#' 
#' Some of the solutions used to build this function were found in the source 
#' code of the R-package \strong{intamapInteractive}. As such, the authors of 
#' that package, Edzer Pebesma <\email{edzer.pebesma@@uni-muenster.de}> and Jon
#' Skoien <\email{jon.skoien@@gmail.com}>, are entitled \sQuote{contributors} to
#' the R-package \pkg{pedometrics}.
#' 
#' @seealso \code{ssaOptim}, \code{\link[sp]{zerodist}}, 
#' \code{\link[base]{jitter}}, \code{\link[geoR]{jitter2d}},
#' \code{\link[rgeos]{gContains}}.
#' @keywords iteration spatial
#' @concept jitter perturb
#' @export
#' @examples
#' require(sp)
#' data(meuse.grid)
#' meuse.grid <- as.matrix(meuse.grid[, 1:2])
#' meuse.grid <- matrix(cbind(1:dim(meuse.grid)[1], meuse.grid), ncol = 3)
#' pts1 <- sample(c(1:dim(meuse.grid)[1]), 155)
#' pts2 <- meuse.grid[pts1, ]
#' pts3 <- spJitterFinite(points = pts2, candidates = meuse.grid, x.min = 40,
#'                       x.max = 100, y.min = 40, y.max = 100, which.point = 10)
#' plot(meuse.grid[, 2:3], asp = 1, pch = 15, col = "gray")
#' points(pts2[, 2:3], col = "red", cex = 0.5)
#' points(pts3[, 2:3], pch = 19, col = "blue", cex = 0.5)
# FUNCTION #####################################################################
spJitterFinite <-
  function (points, candidates, x.max, x.min, y.max, y.min, which.point) {
    pt1 <- .spJitterCpp(points[, 2:3], candidates[, 2:3], x.max, x.min, y.max, 
                        y.min, which.point)
    
    # ASR: Pass all the following to C++
    pt1 <- pt1[pt1 != 0]
    pt2 <- candidates[sample(pt1, 1), ]
    dup <- duplicated(rbind(pt2, points))
    if (any(dup)) {
      while (any(dup)) {
        pt2 <- candidates[sample(pt1, 1), ]
        dup <- duplicated(rbind(pt2, points))
      }
    }
    res <- points
    res[which.point, ] <- pt2
    return (res)
  }
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
# Details
# The second method can be much slower than the first depending on the number of points, on the shape of the area and on how the other arguments are set. This method does not use a finite set of candidate locations. Instead, the number of candidate locations is infinite. Its domain can be defined using the argument \code{where}. The reason for the larger amount of time demanded is that the method has two internal steps to 1) check if the perturbed point falls inside the spatial domain, and b) check if two of more points have coincident coordinates (set using argument \code{zero}). Using an infinite set of candidate locations will usually allow obtaining better results in optimization exercises such as spatial simulated annealing. However, the amount of time may be prohibitive depending on the complexity of the problem.
# 
# The sub-argument \code{max} in both arguments \code{x.coord} and \code{y.coord} defines the lower and upper limits of a uniform distribution (\code{runif(n, min = -max, max = max)}). The quantity of noise added to the coordinates of the point being perturbed is sampled from this uniform distribution. By default, the maximum quantity of random noise added to the x and y coordinates is, respectively, equal to half the width and height of the bounding box of the set of points. This is equivalent to a vector \strong{h} of length equal to half the diagonal of the bounding box. Therefore, a larger jittering is allowed in the longer coordinate axis (x or y).
# 
# The direction of the perturbation is defined by the sign of the values sampled from the uniform distribution. This means that the perturbation can assume any direction from 0 to 360 degrees. By contrast, the function \code{\link[geoR]{jitter2d}} in the R-package \pkg{geoR} samples from a uniform distribution a value for the length of the vector \strong{h} and a value for the direction of the perturbation.
# 
# \code{spJitter} allows to set the minimum quantity of random noise added to a coordinate with the sub-argument \code{min}. The absolute difference between the original coordinate value and the jittered coordinate value is used to evaluate this constraint. If the constraint is not met, \code{min} receives the sign of the value sample from the uniform distribution and is added to the original coordinate value. This does not guarantee that the perturbation will be in the same direction, but in the same quadrant.
# 
# When a spatial domain is defined, \code{spJitter} evaluates if the perturbed points fall inside it using the function \code{\link[rgeos]{gContains}} from the R-package \pkg{rgeos}. All points falling outside the spatial domain are identified and have their original coordinates jittered one again. Every new coordinate falling inside the spatial domain is accepted. Every point falling outside the spatial domain has its coordinates jittered till it falls inside the spatial domain. The number of iterations necessary to meet this constraint depends on the complexity of the shape of the spatial domain. \code{spJitter} tries to speed up the process by linearly decreasing the maximum quantity of noise added to the coordinates at each iteration. If the number of iterations was not enough to guarantee all points inside the spatial domain, \code{spJitter} returns the jittered SpatialPoints with a warning message informing how many points do not meet the constraint.
# spJitterFinite <-
#   function (points, candidates, x.max, x.min, y.max, y.min, which.pts) {
#     d_x <- x.max + x.min
#     d_y <- y.max + y.min
#     pt0 <- points[which.pts, ]
#     d_x <- unlist(c(pt0[1] - d_x, pt0[1] + d_x))
#     d_y <- unlist(c(pt0[2] - d_y, pt0[2] + d_y))
#     pt1 <- which(candidates[, 1] >= d_x[1] & candidates[, 1] <= d_x[2] &
#                    candidates[, 2] >= d_y[1] & candidates[, 2] <= d_y[2])
#     pt2 <- candidates[sample(pt1, 1), ]
#     dup <- duplicated(rbind(pt2, points))
#     if (any(dup)) {
#       while (any(dup)) {
#         pt2 <- candidates[sample(pt1, 1), ]
#         dup <- duplicated(rbind(pt2, points))
#       }
#     }
#     res <- points
#     res[which.pts, ] <- pt2
#     return (res)
#   }
# OLD spJitter FUNCTION ########################################################
# spJitter <- 
#   function (obj, candidates = NULL, where = NULL, which = NULL, finite = NULL,
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
#       if (is.null(candidates)) {
#         stop ("'candidates' is a mandatory argument")
#       }
#       if (!inherits(candidates, what = "data.frame")) {
#         stop ("'candidates' should be a data.frame")
#       } else {
#         if (dim(candidates)[2] != 2 || 
#               any(colnames(candidates) != c("x", "y"))) {
#           stop ("'candidates' should have 2 columns named 'x' and 'y'")
#         }
#       }
#       if (unique(which == "all")) {
#         stop ("this option is not functional yet")
#         #res <- candidates[sample(c(1:length(candidates)), length(obj)), ]
#         } else {
#           if (length(which) > 1) {
#             stop ("this option is not functional yet")
#             #pt0 <- coordinates(obj[which, ])
#             #res <- list()
#             #for (i in 1:length(which)) {
#             #  cand <- coordinates(candidates)
#             #  d_x <- x.coord$max + x.coord$min
#             #  d_y <- y.coord$max + y.coord$min
#             #  d_x <- c(pt0[i, "x"] - d_x, pt0[i, "x"] + d_x)
#             #  d_y <- c(pt0[i, "y"] - d_y, pt0[i, "y"] + d_y)
#             #  cand <- which(cand[, "x"] >= d_x[1] & cand[, "x"] <= d_x[2] &
#             #                  cand[, "y"] >= d_y[1] & cand[, "y"] <= d_y[2])
#             #  res[i] <- candidates[sample(cand, 1), ]
#             #}
#             #res <- data.frame(t(sapply(res, coordinates)))
#             #colnames(res) <- c("x", "y")
#             #coordinates(res) <- ~ x + y
#             #proj4string(res) <- proj4string(obj)
#             #res <- rbind(obj[-which, ], res)
#           } else {
#             res <- .spJitterFiniteOne(obj = obj, candidates = candidates,
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
