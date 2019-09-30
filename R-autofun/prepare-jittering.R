# Prepare for jittering
#
# COMMAND
# eval(.prepare_jittering())
# 
# SUMMARY
# 1. If needed, estimate 'x.max' and 'y.max' from 'candi';
# 2. Set 'x_max0' and 'y_max0';
#
# NOTES
# 1. Estimating 'x.max' and 'y.max' from 'candi' is very expensive. An alternative should be seek.
# 2. A function should be developed so that 'SpatialTools::dist1()' would no longer be a dependency of 
#   'spsann'.

# 'x.min' and 'y.min' are now set to zero by default
x.min <- schedule$x.min
y.min <- schedule$y.min

# Argument `cellsize` of function `scheduleSPSANN` is not estimated from `candi` anymore. This is because 
# `candi` can also be an existing irregular sample configuration and unexperienced users can be unaware of
# the need to set `cellsize = 0` in this case. Only `x.max` and `y.max` are still estimated from `candi`. The
# user now has to inform `cellsize` manually and this should speed up computations.
cellsize <- schedule$cellsize

is_null_x_max <- is.null(schedule$x.max)
is_null_y_max <- is.null(schedule$y.max)
# is_null_cellsize <- is.null(schedule$cellsize)
# if (any(c(is_null_x_max, is_null_y_max, is_null_cellsize) == TRUE)) {
if (any(c(is_null_x_max, is_null_y_max) == TRUE)) {
  # if (!missing(eval.grid)) {
    # message("estimating jittering parameters from 'eval.grid'...")
    # x <- SpatialTools::dist1(as.matrix(eval.grid[, "x"]))
    # y <- SpatialTools::dist1(as.matrix(eval.grid[, "y"]))
  # } else {
    # message("estimating jittering parameters from 'candi'...")
    message("estimating 'x.max' and 'y.max' from 'candi'...")
    x <- SpatialTools::dist1(as.matrix(candi[, "x"]))
    y <- SpatialTools::dist1(as.matrix(candi[, "y"]))
  # }
  x.max <- ifelse(is_null_x_max, max(x) / 2, schedule$x.max)
  # cellsize <- ifelse(is_null_cellsize, min(x[x > 0]), schedule$cellsize)
  y.max <- ifelse(is_null_y_max, max(y) / 2, schedule$y.max)
  # if (is_null_cellsize) { 
    # cellsize <- c(cellsize, min(y[y > 0]))
  # }
  
} else {
  
  # If nothing is missing...
  x.max <- schedule$x.max
  y.max <- schedule$y.max
  # cellsize <- schedule$cellsize
}

x_max0 <- x.max
y_max0 <- y.max

# If 'cellsize' = 0, that is, a FINITE set of candidate locations is used, then we have to compute 'x_min0'
# and 'y_min0'. These are the maximum distance to the nearest neighbouring candidate location in the x- and 
# y-coordinates. Such values are needed so that at the end of the optimization the maximum distance that a 
# given sample point can be jittered is larger than zero. Otherwise it will be stuck in its own location and
# no optimization will happen at all.
# If 'cellsize' > 0, then 'x_min0' and 'y_min0' are set to zero.
# These values are used at the end of each chain when control parameters are updated.
if (all(cellsize) == 0) {
  if (!exists(x = 'x')) {
    x <- SpatialTools::dist1(as.matrix(candi[, "x"]))
    y <- SpatialTools::dist1(as.matrix(candi[, "y"]))
  }
  diag(x) <- Inf
  diag(y) <- Inf
  x_min0 <- max(pedometrics::rowMinCpp(x))
  y_min0 <- max(pedometrics::rowMinCpp(y))
} else {
  x_min0 <- 0
  y_min0 <- 0
}
