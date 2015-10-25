# Set plotting options
#
# SUMMARY
# 1. Check is plotting is required;
# 2. If plotting is required, record the current plotting options;
# 3. Force the plotting options to be reset at the end of the optimization;
# 4. Estimate the 'boundary' if it is missing;
# 5. Open two new plotting devices.
#  
# NOTES
# 1. The estimated boundary is an object of class SpatialPoints. An object of
#    class SpatialPolygon cannot be created because the SpatialPoints are not
#    in the correct order. The following commands could be used otherwise:
#    boundary <- Polygons(list(Polygon(boundary)), ID = as.character(1))
#    boundary <- SpatialPolygons(list(boundary))
# 2. A precise alternative is to use the package 'rgeos' and the following 
#    commands:
#    coordinates(candi) <- ~ x + y
#    gridded(candi) <- TRUE
#    boundary <- as(candi, "SpatialPolygons")
#    boundary <- gUnionCascaded(boundary)
#    The drawback is that the 'rgeos' package would then be a dependency of the
#    'spsann' package, which some users may not like.
#    
if (plotit) {
  par0 <- par()
  on.exit(suppressWarnings(par(par0)))
  if (missing(boundary)) {
    x <- by(candi[, 1], as.factor(candi[, 2]), FUN = range, simplify = FALSE)
    x <- do.call(rbind, x)
    d <- dist(x)
    d <- min(d[d > 0]) / 2
    x[, 1] <- x[, 1] - d
    x[, 2] <- x[, 2] + d
    y <- as.numeric(rownames(x))
    xy <- cbind(c(x[, 1], x[, 2]), rep(y, 2))
    
    y <- by(candi[, 2], as.factor(candi[, 1]), FUN = range, simplify = FALSE)
    y <- do.call(rbind, y)
    d <- dist(y)
    d <- min(d[d > 0]) / 2
    y[, 1] <- y[, 1] - d
    y[, 2] <- y[, 2] + d
    x <- as.numeric(rownames(y))
    yx <- cbind(rep(x, 2), c(y[, 1], y[, 2]))
    
    boundary <- unique(rbind(xy, yx))
    rownames(boundary) <- 1:nrow(boundary)
    boundary <- sp::SpatialPoints(boundary)
  }
  rm(x, d, y, xy, yx)
  
  # Open two new plotting devices
  grDevices::dev.new()
  grDevices::dev.new()
}
#
# COMMAND
# # Set plotting options
# eval(.plotting_options())
