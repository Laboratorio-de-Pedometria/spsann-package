# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
source('R/optimMKV.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')

# 0) DEFAULT EXAMPLE ###########################################################
require(pedometrics)
require(sp)
require(rgeos)
require(gstat)
require(plyr)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(1:dim(candi)[1], candi), ncol = 3)
covars <- as.data.frame(meuse.grid)
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
model <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
set.seed(2001)
res <- optimMKV(points = 100, candi = candi, covars = covars, 
                equation = z ~ dist, model = model, krige.stat = "mean", 
                x.max = x.max, x.min = 40, y.max = y.max, y.min = 40,
                boundary = boundary, iterations = 1000, plotit = TRUE)
tail(attr(res, "energy"), 1) # 11.53921
objMKV(points = res, candi = candi, covars = covars, equation = z ~ dist, 
       model = model, krige.stat = "mean")

# 1) GREEDY ALGORITHM ##########################################################
rm(list = ls())
gc()
require(rgeos)
source('R/optimMKV.R')
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
covars <- as.data.frame(meuse.grid)
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
model <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
set.seed(2001)
res <- optimMKV(points = 100, candi = candi, covars = covars, 
                equation = z ~ dist, model = model, krige.stat = "mean", 
                x.max = x.max, x.min = 40, y.max = y.max, y.min = 40,
                boundary = boundary, iterations = 1000, plotit = TRUE,
                greedy = TRUE)
tail(attr(res, "energy"), 1) # 11.55729
objMKV(points = res, candi = candi, covars = covars, equation = z ~ dist, 
       model = model, krige.stat = "mean")

# 2) MANY COVARIATES ###########################################################
rm(list = ls())
gc()
require(rgeos)
source('R/optimMKV.R')
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
covars <- as.data.frame(meuse.grid)
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
model <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
set.seed(2001)
res <- optimMKV(points = 100, candi = candi, covars = covars, 
                equation = z ~ dist + soil + ffreq + x + y, model = model,
                krige.stat = "mean", x.max = x.max, x.min = 40, y.max = y.max, 
                y.min = 40, boundary = boundary, iterations = 1000, 
                plotit = TRUE)
tail(attr(res, "energy"), 1) # 11.55729
objMKV(points = res, candi = candi, covars = covars, equation = z ~ dist, 
       model = model, krige.stat = "mean")
