# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
source('R/optimMUKV.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
# 0) DEFAULT EXAMPLE ###########################################################
require(pedometrics)
require(sp)
require(rgeos)
require(gstat)
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
y.min <- 40
x.min <- 40
points <- 100
equation <- z ~ dist
model <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
krige.stat <- "mean"
set.seed(2001)
res <- optimMUKV(points = points, candi = candi, covars = covars, 
                 equation = equation, model = model, krige.stat = krige.stat, 
                 x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min,
                 boundary = boundary, iterations = 100, plotit = TRUE)
tail(attr(res, "energy"), 1) # 11.79797
objMUKV(points = res, candi = candi, covars = covars, equation = equation, 
        model = model, krige.stat = krige.stat) # 11.87138
