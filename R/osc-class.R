#' Class \sQuote{OptimizedSampleConfiguration}
#' 
#' An S4 class
#' 
#' @slot points
#' 
#' @slot spsann
#' 
#' @slot objective
#' 
#' @rdname osc-class
#' @export
#' @examples 
#' methods::getSlots("OptimizedSampleConfiguration")
#' methods::is(osc, "OptimizedSampleConfiguration")
# MAIN FUNCTION - OptimizedSampleConfiguration CLASS ###########################
OptimizedSampleConfiguration <-
  methods::setClass(Class = "OptimizedSampleConfiguration", 
                    slots = c(points = "data.frame", 
                              spsann = "list",
                              objective = "list"),
                    
                    prototype = list(
                      
                      # The optimized sample configuration
                      points = data.frame(id = NA_integer_, x = NA_real_, 
                                          y = NA_real_),
                      
                      # Information about the spatial simulated annealing
                      spsann = list(
                        acceptance = data.frame(initial = NA_real_),
                        cellsize = data.frame(x = NA_real_, y = NA_real_),
                        chains = data.frame(total = NA_integer_, 
                                            used = NA_integer_,
                                            length = NA_integer_),
                        jitter = data.frame(x = rep(NA_real_, 2), 
                                            y = rep(NA_real_, 2), 
                                            row.names = c("min", "max")),
                        running = data.frame(time = NA_real_, 
                                             units = NA_character_),
                        stopping = NA_integer_,
                        temperature = data.frame(initial = NA_real_, 
                                                 final = NA_real_)),
                      
                      # Information about the objective function
                      objective = list(
                        name = NA_character_,
                        energy = data.frame(NA_real_),
                        nadir = data.frame(NA_real_),
                        utopia = data.frame(NA_real_),
                        weights = data.frame(NA_real_))
                    ))
