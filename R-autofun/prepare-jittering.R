# Prepare for jittering
#
# COMMAND
# eval(.prepare_jittering())
# 
# SUMMARY
# 1. Check if 'x.min', 'x.max', 'y.min', and 'y.max', are missing;
# 2. If missing, estimate them from 'candi';
# 3. Set 'x_max0' and 'y_max0', which are used later during the optimization;
# 4. FUTURE FEATURE: Determine 'cellsize', the x- and y- dimensions of the
#    individual elements of 'candi'.
#
# NOTES
# 1. Estimating 'x.min', 'x.max', 'y.min', and 'y.max' from 'candi' is very
#    expensive. An alternative should be seek.
# 2. A function should be developed so that 'SpatialTools::dist1()' would no
#    longer be a dependency of 'spsann'.
# 3. 'cellsize' will be used in the future when the use of an infinite set of
#    candidate locations will be enabled. The idea is to use the same strategy
#    used in 'spcosa::spsample()'.
#
if (missing(x.min) && missing(x.max) && missing(y.min) && missing(y.max)) {
  message("estimating 'x.min', 'x.max', 'y.min', and 'y.max' from 'candi'")
  x <- SpatialTools::dist1(as.matrix(candi[, "x"]))
  id <- x > 0
  x.min <- min(x[id])
  x.max <- max(x) / 2
  
  y <- SpatialTools::dist1(as.matrix(candi[, "y"]))
  id <- y > 0
  y.min <- min(y[id])
  y.max <- max(y) / 2
  
  rm(x, id, y)
}
x_max0 <- x.max
y_max0 <- y.max
# if (missing(cellsize)) {
#   cellsize <- c(x.min, y.min)
# }
