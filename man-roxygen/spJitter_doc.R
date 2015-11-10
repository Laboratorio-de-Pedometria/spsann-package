#  Template documentation for spatial jittering
################################################################################
#' @section Jittering methods:
#' There are two ways of jittering the coordinates. They differ on how the
#' set of candidate locations is defined. The first method uses an 
#' \emph{infinite} set of candidate locations, that is, any locations in the 
#' spatial domain can be selected as the new location of a jittered point. All 
#' that this method needs is a polygon indicating the boundary of the spatial 
#' domain. This method is not implemented in the \pkg{spsann}-package because 
#' it is computationally demanding: every time a point is jittered, it is 
#' necessary to check if the point is inside the spatial domain.
#' 
#' The second method consists of using a \emph{finite} set of candidate 
#' locations for the jittered points. A finite set of candidate locations is
#' created by discretizing the spatial domain, that is, creating a fine grid of
#' points that serve as candidate locations for the jittered points. This is
#' the least computationally demanding jittering method because, by definition, 
#' the jittered point will always be inside the spatial domain. However, not 
#' all locations in the spatial domain can be selected as the new location for 
#' a jittered point.
#' 
#' Using a finite set of candidate locations has another important 
#' inconvenience. When a point is jittered, it may be that the new location 
#' already is occupied by another point. If this happens, another location has
#' to be iteratively sought for, say, as many times as the number of points in
#' the system. In general, the more points there are in the system, the more 
#' likely it is that the new location already is occupied by another point. If 
#' a solution is not found in a reasonable time, the point selected to be 
#' jittered is kept in its original location. Such a proceedure is clearly 
#' suboptimal.
#' 
#' The \pkg{spsann}-package uses a more elegant method which is based on using 
#' a finite set of candidate locations coupled with a form of \emph{two-stage 
#' random sampling} as implemented in \code{\link[spcosa]{spsample}}. Because 
#' the candidate locations are placed on a finite regular grid, they can be 
#' seen as being the centre nodes of a finite set of grid cells (or pixels of 
#' a raster image). In the first stage, one of the \dQuote{grid cells} is 
#' selected with replacement, i.e. independently of already being occupied by 
#' another sample point. The new location for the point chosen to be jittered 
#' is selected within that \dQuote{grid cell} by simple random sampling. This 
#' method guarantees that \emph{almost} any location in the spatial domain can 
#' be a candidate location. It also discards the need to check if the new 
#' location already is occupied by another point, simplifying the code and 
#' (possibly) speeding up the computations.
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

