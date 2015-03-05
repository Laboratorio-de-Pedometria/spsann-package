# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
source('R/optimCORR.R')
source('R/optimACDC.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
# 0) DEFAULT EXAMPLE ###########################################################
require(pedometrics)
require(sp)
require(rgeos)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(1:dim(candi)[1], candi), ncol = 3)
covars <- meuse.grid[, 5]
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
y.min <- 40
x.min <- 40
iterations <- 1000
points <- 100
use.coords <- TRUE
set.seed(2001)
res <- optimCORR(points = points, candi = candi, covars = covars, 
                 use.coords = use.coords, x.max = x.max, x.min = x.min, 
                 y.max = y.max, y.min = y.min, boundary = boundary, 
                 iterations = iterations)
tail(attr(res, "energy"), 1) # 0.001747814
objCORR(points = res, candi = candi, covars = covars, use.coords = use.coords)
# 1) FACTOR COVARIATES WITH THE COORDINATES ####################################
rm(list = ls())
gc()
source('R/optimCORR.R')
source('R/optimACDC.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
data(meuse.grid)
candi <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(1:dim(candi)[1], candi), ncol = 3)
covars <- meuse.grid[, 6:7]
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
y.min <- 40
x.min <- 40
strata.type <- "area"
use.coords <- TRUE
iterations <- 1000
points <- 100
set.seed(2001)
res <- optimCORR(points = points, candi = candi, covars = covars, 
                 strata.type = strata.type, use.coords = use.coords,
                 x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, 
                 boundary = boundary, iterations = iterations)
tail(attr(res, "energy"), 1) # 1.249576
objCORR(points = res, candi = candi, covars = covars, use.coords = use.coords)
# 2) FACTOR COVARIATES #########################################################
rm(list = ls())
gc()
source('R/optimCORR.R')
source('R/optimACDC.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
data(meuse.grid)
candi <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(1:dim(candi)[1], candi), ncol = 3)
str(meuse.grid)
covars <- meuse.grid[, rep(6:7, 3)]
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
y.min <- 40
x.min <- 40
use.coords <- FALSE
iterations <- 1000
points <- 100
set.seed(2001)
res <- optimCORR(points = points, candi = candi, covars = covars, 
                 use.coords = use.coords, x.max = x.max, x.min = x.min, 
                 y.max = y.max, y.min = y.min, 
                 boundary = boundary, iterations = iterations)
tail(attr(res, "energy"), 1) # 0.01085086
objCORR(points = res, candi = candi, covars = covars, use.coords = use.coords)
