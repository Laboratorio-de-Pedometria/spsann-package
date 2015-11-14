OptimizedSampleConfiguration <-
  setClass(Class = "OptimizedSampleConfiguration", 
           slots = c(points = "data.frame", 
                     spsann = "list",
                     objective = "list"),
           
           prototype = list(
             
             # The optimized sample configuration
             points = data.frame(id = NA_integer_, x = NA_real_, y = NA_real_),
             
             # Information about the spatial simulated annealing
             spsann = list(
               acceptance = data.frame(initial = NA_real_),
               cellsize = data.frame(x = NA_real_, y = NA_real_),
               chains = data.frame(total = NA_integer_, used = NA_integer_,
                                   length = NA_integer_),
               jitter = data.frame(x = rep(NA_real_, 2), y = rep(NA_real_, 2), 
                                   row.names = c("min", "max")),
               running = data.frame(time = NA_real_, units = NA_character_),
               stopping = NA_integer_,
               temperature = data.frame(initial = NA_real_, final = NA_real_)),
             
             # Information about the objective function
             objective = list(
               energy = data.frame(NA_real_),
               name = NA_character_)
             ))
# tmp <- new("OptimizedSampleConfiguration")
# slot(tmp, "objective") <- list(energy = data.frame(NA_real_), name = "DIST")
# str(tmp, 3)
