# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
source('R/optimDIST.R')
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
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, x.max = x.max, x.min = 40, y.max = y.max, 
                 y.min = 40, boundary = boundary, iterations = 100)
tail(attr(res, "energy"), 1) # 1.656926
objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 1) GREEDY ALGORITHM ##########################################################
rm(list = ls())
gc()
source('R/optimDIST.R')
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
covars <- meuse.grid[, 5]
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, greedy = TRUE,
                 use.coords = TRUE, x.max = x.max, x.min = 40, y.max = y.max, 
                 y.min = 40, boundary = boundary, iterations = 100)
tail(attr(res, "energy"), 1) # 1.530403
objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 2) FACTOR COVARIATES WITH THE COORDINATES ####################################
rm(list = ls())
gc()
source('R/optimDIST.R')
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
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
covars <- meuse.grid[, 6:7]
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, 
                 strata.type = "area", use.coords = TRUE, x.max = x.max, 
                 x.min = 40, y.max = y.max, y.min = 40, boundary = boundary,
                 iterations = 100)
tail(attr(res, "energy"), 1) # 1.188592
objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE,
        strata.type = "area")

# 3) FACTOR AND NUMERIC COVARIATES WITH THE COORDINATES ########################
rm(list = ls())
gc()
source('R/optimDIST.R')
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
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
covars <- meuse.grid[, c(1, 2, 5:7)]
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, 
                 strata.type = "area", use.coords = TRUE, x.max = x.max, 
                 x.min = 40, y.max = y.max, y.min = 40, boundary = boundary, 
                 iterations = 100)
tail(attr(res, "energy"), 1) # 2.958176
objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE, 
        strata.type = "area")
