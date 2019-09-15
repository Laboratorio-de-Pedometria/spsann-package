# Prepare for jittering
#
# COMMAND
# eval(.prepare_jittering())
# 
# SUMMARY
# 1. If needed, estimate 'x.min', 'x.max', 'y.min', 'y.max' and 'cellseze' from 'candi' or 'eval.grid';
# 2. Set 'x_max0' and 'y_max0';
#
# NOTES
# 1. Estimating 'x.min', 'x.max', 'y.min', 'y.max' and 'cellsize' from 'candi' or 'eval.grid' is very 
#    expensive. An alternative should be seek.
# 2. A function should be developed so that 'SpatialTools::dist1()' would no longer be a dependency of 
#   'spsann'.

x.min <- schedule$x.min
y.min <- schedule$y.min
is_null_x_max <- is.null(schedule$x.max)
is_null_y_max <- is.null(schedule$y.max)
is_null_cellsize <- is.null(schedule$cellsize)
if (any(c(is_null_x_max, is_null_y_max, is_null_cellsize) == TRUE)) {
  if (!missing(eval.grid)) {
    message("estimating jittering parameters from 'eval.grid'...")
    x <- SpatialTools::dist1(as.matrix(eval.grid[, "x"]))
    y <- SpatialTools::dist1(as.matrix(eval.grid[, "y"]))
  } else {
    message("estimating jittering parameters from 'candi'...")
    x <- SpatialTools::dist1(as.matrix(candi[, "x"]))
    y <- SpatialTools::dist1(as.matrix(candi[, "y"]))
  }
  x.max <- ifelse(is_null_x_max, max(x) / 2, schedule$x.max)
  cellsize <- ifelse(is_null_cellsize, min(x[x > 0]), schedule$cellsize)
  y.max <- ifelse(is_null_y_max, max(y) / 2, schedule$y.max)
  if (is_null_cellsize) { 
    cellsize <- c(cellsize, min(y[y > 0]))
  }
  
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
