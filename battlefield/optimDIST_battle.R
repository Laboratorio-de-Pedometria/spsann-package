# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
source('R/optimDIST.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
source("/home/alessandro/PROJECTS/r-packages/pedometrics/R/numint.R")
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
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, x.max = x.max, x.min = x.min, 
                 y.max = y.max, y.min = y.min, boundary = boundary, 
                 iterations = 100)
tail(attr(res, "energy"), 1) # 165.3155
objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)
# 1) FACTOR COVARIATES WITH THE COORDINATES ####################################
rm(list = ls())
gc()
source("/home/alessandro/PROJECTS/r-packages/pedometrics/R/numint.R")
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
source('R/optimDIST.R')
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
y.min <- 40
x.min <- 40
strata.type <- "area"
covars <- meuse.grid[, 6:7]
use.coords <- TRUE
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, 
                 strata.type = "area", use.coords = use.coords, x.max = x.max, 
                 x.min = x.min, y.max = y.max, y.min = y.min, 
                 boundary = boundary, iterations = 100)
tail(attr(res, "energy"), 1) # 125.0564
objDIST(points = res, candi = candi, covars = covars, use.coords = use.coords,
        strata.type = strata.type)
# 2) FACTOR AND NUMERIC COVARIATES WITH THE COORDINATES ########################
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
source('R/optimDIST.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
source("/home/alessandro/PROJECTS/r-packages/pedometrics/R/numint.R")
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
y.min <- 40
x.min <- 40
points <- 100
iterations <- 100
strata.type <- "area"
covars <- meuse.grid[, c(1, 2, 5:7)]
str(covars)
use.coords <- TRUE
set.seed(2001)
res <- optimDIST(points = points, candi = candi, covars = covars, strata.type = strata.type, use.coords = use.coords, x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, boundary = boundary, iterations = iterations)
tail(attr(res, "energy"), 1) # 305.554
objDIST(points = res, candi = candi, covars = covars, use.coords = use.coords, 
        strata.type = strata.type)
