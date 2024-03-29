#  Template documentation for the cooling schedule
####################################################################################################
#' @details
#' \subsection{Annealing schedule}{
#' The _search graph_ corresponds to the set of effective candidate locations for a sample location
#' selected to be jittered. The size of the search graph, i.e. area within which a sample location
#' can be moved around, is related to the concept of _temperature_. A larger search graph is
#' equivalent to higher temperatures, which potentially result in more movement -- or
#' \sQuote{agitation} -- of the set of sample locations.
#'
#' The current version of the \pkg{spsann}-package uses a linear cooling schedule which depends upon
#' the number of jitters to control the size of the search graph. The equations are
#'
#' \eqn{x_max = x_max0 - (chains_i / chains) * (x_max0 - x_min) + x_cellsize + x_min0}
#'
#' and
#'
#' \eqn{y_max = y_max0 - (chains_i / chains) * (y_max0 - y_min) + y_cellsize + y_min0},
#'
#' where $x_max0$ and $y_max0$ are the maximum allowed shifts in the x- and y-coordinates in the
#' first chain, $x_min$ and $y_min$ are the minimum required shifts in the x- and y-coordinates,
#' $x_max$ and $y_max$ are the maximum allowed shifts in the x- and y-coordinates during the next
#' chain, $chains$ and $chain_i$ are the total and current chains, and $x_cellsize$ and $y_cellsize$
#' are the grid spacing in the x- and y-coordinates. Because $x_cellsize$ and $y_cellsize$ can be
#' equal to zero when a finite set of candidate locations is used, $x_min0$ and $y_min0$ are the
#' maximum nearest neighbour distance in the x- and y-coordinates between candidate locations.
#' }
