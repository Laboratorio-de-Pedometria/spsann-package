#  Template documentation for spatial jittering (include in functions other than spJitter)
###############################################################################################################
#' @details
#' \subsection{Generating mechanism}{
#' There are multiple mechanism to generate a new sample configuration out of an existing one. The main step
#' consists of randomly perturbing the coordinates of a single sample, a process known as \sQuote{jittering}.
#' These mechanisms can be classified based on how the set of candidate locations for the samples is defined.
#' For example, one could use an _infinite_ set of candidate locations, that is, any location in the spatial
#' domain can be selected as a new sample location after a sample is jittered. All that is needed is a 
#' polygon indicating the boundary of the spatial domain. This method is more computationally demanding 
#' because every time an existing sample is jittered, it is necessary to check if the new sample location 
#' falls in spatial domain.
#'
#' Another approach consists of using a _finite_ set of candidate locations for the samples. A finite set of
#' candidate locations is created by discretising the spatial domain, that is, creating a fine (regular) grid
#' of points that serve as candidate locations for the jittered sample. This is a less computationally
#' demanding jittering method because, by definition, the new sample location will always fall in the spatial
#' domain.
#' 
#' Using a finite set of candidate locations has two important inconveniences. First, not all locations in 
#' the spatial domain can be selected as the new location for a jittered sample. Second, when a sample is 
#' jittered, it may be that the new location already is occupied by another sample. If this happens, another 
#' location has to be iteratively sought for, say, as many times as the size of the sample configuration. In 
#' general, the larger the size of the sample configuration, the more likely it is that the new location 
#' already is occupied by another sample. If a solution is not found in a reasonable time, the the sample 
#' selected to be jittered is kept in its original location. Such a procedure clearly is suboptimal.
#' 
#' __spsann__ uses a more elegant method which is based on using a finite set of candidate locations coupled
#' with a form of _two-stage random sampling_ as implemented in the R-package
#' \code{[spsample](https://CRAN.R-project.org/package=spcosa)}. Because the candidate locations are placed on
#' a finite regular grid, they can be taken as the centre nodes of a finite set of grid cells (or pixels of a 
#' raster image). In the first stage, one of the \dQuote{grid cells} is selected with replacement, i.e. 
#' independently of already being occupied by another sample. The new location for the sample chosen to be
#' jittered is selected within that \dQuote{grid cell} by simple random sampling. This method guarantees that 
#' virtually any location in the spatial can be selected. It also discards the need to check if the new 
#' location already is occupied by another sample, speeding up the computations when compared to the first two
#' approaches.
#' }
#'  
#' @note
#' \subsection{Distance between two points}{
#' __spsann__ always computes the distance between two locations (points) as the 
#' [Euclidean distance](https://en.wikipedia.org/wiki/Euclidean_distance) between them. This computation 
#' requires the optimization to operate in the two-dimensional Euclidean space, i.e. the coordinates of
#' the sample, candidate and evaluation locations must be planar coordinates, generally in metres or 
#' kilometres. __spsann__ has no mechanism to check if the coordinates are planar: you are the sole responsible
#' for making sure that this requirement is attained.
#' }
