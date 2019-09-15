#  Template documentation for the cooling schedule
###############################################################################################################
#' @details
#' \subsection{Annealing schedule}{
#' The _search graph_ corresponds to the set of effective candidate locations for a point selected to be
#' jittered. The size of the search graph, i.e. area within which a point can be moved around, is correlated 
#' with the concept of _temperature_. A larger search graph is equivalent to higher temperatures, which
#' potentially result in more movement or \sQuote{agitation} of the set of points or \sQuote{particles}.
#' 
#' The current version of the \pkg{spsann}-package uses a linear cooling schedule which depends upon the 
#' number of jitters to control the size of the search graph. The equations are
#' 
#' \eqn{x_max = x_max0 - (chains_i / chains) * (x_max0 - x_min) + x_cellsize}
#' 
#' and
#' 
#' \eqn{y_max = y_max0 - (chains_i / chains) * (y_max0 - y_min) + y_cellsize},
#' 
#' where $x_max0$ and $y_max0$ are the maximum allowed shifts in the x- and y-coordinates in the first chain, 
#' $x_min$ and $y_min$ are the minimum required shifts in the x- and y-coordinates, $x_max$ and $y_max$ are 
#' the maximum allowed shifts in the x- and y-coordinates during the next chain, $chains$ and $chain_i$ are 
#' the total and current chains, and $x_cellsize$ and $y_cellsize$ are the grid spacing in the x- and 
#' y-coordinates.
#' }
