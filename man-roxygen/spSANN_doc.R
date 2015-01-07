#  Template documentation for spatial simulated annealing
################################################################################
#' @param iterations Integer value defining the maximum number of iterations 
#' that should be used for the optimization. See \sQuote{Details} for more 
#' information.
#' 
#' @param acceptance List with two sub-arguments: \code{initial} and 
#' \code{cooling}. \code{initial} is a numeric value between 0 and 1 defining
#' the initial acceptance probability. Defaults to \code{initial = 0.99}.
#' \code{cooling} is a numeric value defining the exponential factor by with 
#' the acceptance probability decreases at each iteration. Defaults to 
#' \code{cooling = iterations / 10}. See \sQuote{Details} for more 
#' information.
#' 
#' @param stopping A list with one sub-argument: \code{max.count}. 
#' \code{max.count} is an integer value defining the maximum allowable number 
#' of iterations without improvement of the objective function value. This is 
#' also known as the freezing criterion. Defaults to 
#' \code{max.count = iterations / 10}. See \sQuote{Details} for more 
#' information.
#' 
#' @param plotit Logical value for ploting the optimization results. This 
#' includes a) the progress of the objective function values and acceptance 
#' probabilities, and b) the original points, the perturbed points and the 
#' progress of the maximum perturbation in the x and y coordinates. The plots 
#' are updated at each 10 iterations. The boundary of the spatial domain is 
#' passed using the argument \code{boundary}. Defaults to \code{plotit = TRUE}.
#' 
#' @param boundary SpatialPolygon defining the boundary of the spatial domain. 
#' It is mandatory if \code{plotit = TRUE}.
#' 
#' @param progress Logical value for printing a progress bar. Defaults to 
#' \code{progress = TRUE}.
#' 
#' @param verbose Logical value for printing messages about the progress of the
#' optimization.
#' 
#' @section Spatial simulated annealing:
#' \subsection{Search graph}{
#' The search graph corresponds to the set of effective candidate locations for
#' a point being jittered in a given iteration. The size of the search graph, 
#' i.e. the maximum distance that a point can be moved around, is correlated 
#' with the concept of \strong{temperature}. A larger search graph is equivalent
#' to higher temperatures, which potentially result in more movement or 
#' \dQuote{agitation} of the set of points or \dQuote{particles}.
#' 
#' The current implementation of spatial simulated annealing uses a 
#' \strong{linear cooling schedule} depending upon the number of iterations to 
#' control the size of the search graph. The equations are as follows:
#' 
#' \verb{
#' x.max.b <- x.max.a - k / iterations * (x.max.a - x.min)
#' y.max.b <- y.max.a - k / iterations * (y.max.a - y.min)
#' }
#' 
#' where \code{x.max.a} and \code{y.max.a} are the maximum allowed shift in the
#' x and y coordinates in the current iteration, \code{x.min} and \code{y.min} 
#' are the minimum required shift in the x and y coordinates, and \code{x.max.b}
#' and \code{y.max.b} are the maximum allowed shift in the x and y coordinates
#' in the next iteration. \code{iterations} is the total number of iterations 
#' and \code{k} is the current iteration.
#' }
#' \subsection{Acceptance probability}{
#' The acceptance probability is the chance of accepting a new system 
#' configuration that is worse than the current system configuration. The 
#' concept of acceptance probability is related with that of 
#' \strong{temperature}. A higher acceptance probability is equivalent to higher
#' temperatures, which potentially result in more movement or 
#' \dQuote{agitation} of the set of points or \dQuote{particles}.
#' 
#' Using a low initial acceptance probability turns the spatial simulated 
#' annealing into a \emph{greedy} algorithm. It will converge in a shorter time,
#' but the solution found is likely to be a local optimum instead of the global 
#' optimum. Using a high initial acceptance probability (\code{>0.8}) usually is
#' the wisest choice.
#' 
#' An \strong{exponential cooling schedule} depending upon the number of 
#' iterations is used in the current implementation of the spatial simulated
#' annealing to control the acceptance probability. The acceptance probability
#' at each iteration is calculates as follows:
#' 
#' \verb{actual_prob <- acceptance$initial * exp(-k / acceptance$cooling)}
#' 
#' where \code{actual_prob} is the acceptance probability at the \code{k}-th 
#' iteration, \code{acceptance$initial} is the initial acceptance probability, 
#' and \code{acceptance$cooling} is the exponential cooling factor.
#' }
#' \subsection{Starting system configuration}{
#' Unidimensional criterion such as the number of points per lag distance class
#' are dependent on the starting system configuration by definition. This means
#' that, depending on the parameters passed to the spatial simulated annealing 
#' algorithm, many points will likely to stay close to their starting positions.
#' It would be reasonable to use a starting system configuration that is close 
#' to the global optimal, but such thing is not feasible.
#' 
#' Increasing the initial acceptance probability does not guarantee the 
#' independence from the starting system configuration. The most efficient 
#' option in the current implementation of the spatial simulated annealing
#' algorithm is to start using the entire spatial domain as search graph. This 
#' is set using the interval of the x and y coodinates to set \code{x.max} 
#' and \code{y.max} (See above).
#' 
#' An alternative is to start jittering (randomly perturbing) several points at
#' a time and use a cooling schedule to \strong{exponentially} decrease the 
#' number of points jittered at each iteration. The current implementation of 
#' the spatial simulated annealing does not explore such alternative. The 
#' cooling schedule would be as follows:
#' 
#' \verb{
#' new.size <- round(c(old.size - 1) * exp(-k / size.factor) + 1)
#' }
#' 
#' where \code{old.size} and \code{new.size} are the number of points jittered
#' in the previous and next iterations, \code{size.factor} is the cooling 
#' parameter, and \code{k} is the number of the current iteration. The larger 
#' the difference between the starting system configuration and the global 
#' optimum, the larger the number of points that would need to be jittered in
#' the first iterations. This will usually increase the time spent on the first
#' iterations.
#' }
#' \subsection{Number of iterations}{
#' The number of iterations has a large influence on the performance of the 
#' spatial simulated annealing algorithm. The larger the number of possible 
#' system configurations, the higher should the number of iterations be.
#' 
#' The number of possible system configurations increases with:
#' \itemize{
#' \item a high initial acceptance probability
#' \item the use of an infinite set of candidate locations
#' \item the use of a very dense finite set of candidate locations
#' }
#' }
#' @keywords spatial optimize
