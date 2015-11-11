# Prepare for jittering
#
# COMMAND
# eval(.prepare_jittering())
# 
# SUMMARY
# 1. If needed, estimate 'x.min', 'x.max', 'y.min', and 'y.max' from 'candi';
# 2. Set 'x_max0' and 'y_max0';
# 3. FUTURE FEATURE: Determine 'cellsize', the x- and y- dimensions of the
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

x.min <- schedule$x.min
y.min <- schedule$y.min

aa <- is.null(schedule$x.max)
bb <- is.null(schedule$y.max)
cc <- is.null(schedule$cellsize)
if (any(c(aa, bb, cc) == TRUE)) {
  
  message("estimating jittering parameters from 'candi'...")
  x <- SpatialTools::dist1(as.matrix(candi[, "x"]))
  id <- x > 0
  x.max <- ifelse(aa, max(x) / 2, schedule$x.max)
  if (cc) { cellsize <- min(x[id]) } else { cellsize <- schedule$cellsize }
  
  y <- SpatialTools::dist1(as.matrix(candi[, "y"]))
  id <- y > 0
  y.max <- ifelse(bb, max(y) / 2, schedule$y.max)
  if (cc) { cellsize <- c(cellsize, min(y[id])) }
  
} else {
  
  # If nothing is missing...
  x.max <- schedule$x.max
  y.max <- schedule$y.max
  cellsize <- schedule$cellsize
}

x_max0 <- x.max
y_max0 <- y.max
# if (is.null(schedule$x.min) && is.null(schedule$x.max) && 
#     is.null(schedule$y.min) && is.null(schedule$y.max)) {
#   message("estimating 'x.min', 'x.max', 'y.min', and 'y.max' from 'candi'")
#   x <- SpatialTools::dist1(as.matrix(candi[, "x"]))
#   id <- x > 0
#   x.min <- min(x[id])
#   x.max <- max(x) / 2
#   
#   y <- SpatialTools::dist1(as.matrix(candi[, "y"]))
#   id <- y > 0
#   y.min <- min(y[id])
#   y.max <- max(y) / 2
#   
#   rm(x, id, y)
# } else {
#   x.min <- schedule$x.min
#   x.max <- schedule$x.max
#   y.min <- schedule$y.min
#   y.max <- schedule$y.max
# }
# x_max0 <- x.max
# y_max0 <- y.max
# if (is.null(schedule$cellsize)) {
#   schedule$cellsize <- c(x.min, y.min)
# }
