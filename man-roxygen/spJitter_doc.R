#  Template documentation for spatial jittering
################################################################################
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
#' perturbed points. \code{candi} must have two columns in the following 
#' order: \code{[, "x"]} the projected x-coordinates, and \code{[, "y"]} the 
#' projected y-coordinates.
#'
#' @param x.max,x.min,y.max,y.min Numeric value. The minimum and maximum 
#' quantity of random noise to be added to the projected x- and y-coordinates.
#' The minimum quantity should be equal to, at least, the minimum distance
#' between two neighbouring candidate locations. The units are the same as of 
#' the projected x- and y-coordinates. If missing, they are estimated from
#' \code{candi}.
#'
#' @section Jittering methods:
#' There are two ways of jittering the coordinates. They differ on how the
#' set of candidate locations is defined. The first method uses an 
#' \emph{infinite} set of candidate locations, that is, any point in the spatial
#' domain can be selected as the new location of a perturbed point. All that 
#' this method needs is a polygon indicating the boundary of the spatial domain.
#' This method is not implemented in the \pkg{spsann} package (yet) because it 
#' is computationally demanding: every time a point is jittered, it is necessary
#' to check if it is inside the spatial domain.
#' 
#' The second method consists of using a \emph{finite} set of candidate 
#' locations for the perturbed points. A finite set of candidate locations is
#' created by discretizing the spatial domain, that is, creating a fine grid of
#' points that serve as candidate locations for the perturbed points. This is 
#' the only method currently implemented in the \pkg{spsann} package because it
#' is one of the least computationally demanding.
#' 
#' Using a finite set of candidate locations has one important inconvenience.
#' When a point is selected to be jittered, it may be that the new location 
#' already is occupied by another point. If this happens, another location is 
#' iteratively sought for as many times as there are points in \code{points}. 
#' Because the more points there are in \code{points}, the more likely it is 
#' that the new location already is occupied by another point. If a solution is
#' not found, the point selected to be jittered point is kept in its original 
#' location.
#' 
#' A more elegant method can be defined using a finite set of candidate 
#' locations coupled with a form of \emph{two-stage random sampling} as 
#' implemented in \code{\link[spcosa]{spsample}}. Because the candidate 
#' locations are placed on a finite regular grid, they can be seen as being the 
#' centre nodes of a finite set of grid cells (or pixels of a raster image). In 
#' the first stage, one of the \dQuote{grid cells} is selected with replacement,
#' i.e. independently of already being occupied by another sample point. The new
#' location for the point chosen to be jittered is selected within that  
#' \dQuote{grid cell} by simple random sampling. This method guarantees that 
#' any location in the spatial domain can be a candidate location. It also
#' discards the need to check if the new location already is occupied by
#' another point. This method is not implemented (yet) in the \pkg{spsann} 
#' package.
#'
#' @section Distance between two points:
#' The distance between two points is computed as the Euclidean distance between
#' them. This computation assumes that the optimization is operating in the 
#' two-dimensional Euclidean space, i.e. the coordinates of the sample points 
#' and candidate locations should not be provided as latitude/longitude. Package 
#' \pkg{spsann} has no mechanism to check if the coordinates are projected, and
#' the user is responsible for making sure that this requirement is attained.
#'
#' @references
#' Edzer Pebesma, Jon Skoien with contributions from Olivier Baume, A. Chorti, 
#' D.T. Hristopulos, S.J. Melles and G. Spiliopoulos (2013). 
#' \emph{intamapInteractive: procedures for automated interpolation - methods 
#' only to be used interactively, not included in \code{intamap} package.} R 
#' package version 1.1-10.
#' 
#' van Groenigen, J.-W. \emph{Constrained optimization of spatial sampling: 
#' a geostatistical approach.} Wageningen: Wageningen University, p. 148, 1999.

