# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(SpatialTools)
require(sp)
require(rgeos)
source('R/spSANNtools.R')
source('R/spJitter.R')
source('R/optimMSSD.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
Rcpp::sourceCpp('src/calcMSSDCpp.cpp')
Rcpp::sourceCpp('src/updateMSSDCpp.cpp')
# 0) DEFAULT EXAMPLE ###########################################################
require(pedometrics)
require(sp)
require(rgeos)
require(SpatialTools)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(c(1:dim(candi)[1]), candi), ncol = 3)
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
x.min <- 40
y.min <- 40
iterations <- 1000
points <- 100
set.seed(2001)
res <- optimMSSD(points = points, candi = candi, x.max = x.max, x.min = x.min,
                 y.max = y.max, y.min = y.min, iterations = iterations,
                 boundary = boundary)
tail(attr(res, "energy.state"), 1) # 9896.487
objMSSD(candi = candi, points = res)
