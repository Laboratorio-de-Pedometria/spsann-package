# Generated by autofun (0.0.0.9000): do not edit by hand!!!
# Please edit source code in spsann-package/R-autofun/prepare-jittering.R
.prepare_jittering<-function(...){
expression(x.min <- schedule$x.min, y.min <- schedule$y.min, 
    cellsize <- schedule$cellsize, x.max <- ifelse(test = is.null(schedule$x.max), yes = diff(range(candi[, 'x'])) / 2, no = schedule$x.max), 
    y.max <- ifelse(test = is.null(schedule$y.max), yes = diff(range(candi[, 'y'])) / 2, no = schedule$y.max), 
    x_max0 <- x.max, y_max0 <- y.max, if (all(cellsize == 0)) {
      x_min0 <- SpatialTools::dist1(as.matrix(candi[, "x"]))
      y_min0 <- SpatialTools::dist1(as.matrix(candi[, "y"]))
      diag(x_min0) <- Inf
      diag(y_min0) <- Inf
      x_min0 <- max(pedometrics::rowMinCpp(x_min0))
      y_min0 <- max(pedometrics::rowMinCpp(y_min0))
    } else {
      x_min0 <- 0
      y_min0 <- 0
    })
}

