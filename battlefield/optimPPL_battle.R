# Initial settings
rm(list = ls())
gc()
source('R/spSANNtools.R')
source('R/spJitter.R')
source('R/optimPPL.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
Rcpp::sourceCpp('src/updatePPLCpp.cpp')
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
candi <- matrix(cbind(1:nrow(candi), candi), ncol = 3)
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
x.min <- 40
y.min <- 40
cutoff <- sqrt((x.max * x.max) + (y.max * y.max)) / 2
iterations <- 50000
points <- 400
lags <- 7
lags.base <- 2
criterion <- "distribution"
lags.type <- "exponential"
pairs <- FALSE
set.seed(2001)
res <- optimPPL(points = points, candi = candi, lags = lags, pairs = pairs,
                lags.base = lags.base, criterion = criterion, cutoff = cutoff,
                lags.type = lags.type,  x.max = x.max, x.min = x.min, 
                y.max = y.max, y.min = y.min, boundary = boundary,
                iterations = iterations, plotit = FALSE, verbose = FALSE)
countPPL(points = res, lags = lags, lags.type = lags.type, pairs = pairs,
         lags.base = lags.base, cutoff = cutoff)
tail(attr(res, "energy.state"), 1) # 65
objPPL(points = res, lags = lags, lags.type = lags.type, pairs = pairs,
       lags.base = lags.base, cutoff = cutoff, criterion = criterion)
# 1) Point pairs ###############################################################
rm(list = ls())
gc()
source('R/spSANNtools.R')
source('R/spJitter.R')
source('R/optimPPL.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
Rcpp::sourceCpp('src/updatePPLCpp.cpp')
data(meuse.grid)
candi <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(1:nrow(candi), candi), ncol = 3)
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
x.min <- 40
y.min <- 40
cutoff <- sqrt((x.max * x.max) + (y.max * y.max))
iterations <- 100
points <- 100
lags <- 7
lags.base <- 2
criterion <- "distribution"
lags.type <- "exponential"
pairs <- TRUE
set.seed(2001)
res <- optimPPL(points = points, candi = candi, lags = lags, pairs = pairs,
                lags.base = lags.base, criterion = criterion, cutoff = cutoff,
                lags.type = lags.type,  x.max = x.max, x.min = x.min, 
                y.max = y.max, y.min = y.min, boundary = boundary,
                iterations = iterations)
countPPL(points = res, lags = lags, lags.type = lags.type, 
         lags.base = lags.base, cutoff = cutoff, pairs = pairs)
tail(attr(res, "energy.state"), 1) # 7592.857
objPPL(points = res, lags = lags, lags.type = lags.type, pairs = pairs,
       lags.base = lags.base, cutoff = cutoff, criterion = criterion)

# 2) Points per lag - select sample points from candi ##########################
rm(list = ls())
gc()
source('R/spSANNtools.R')
source('R/spJitter.R')
source('R/optimPPL.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
Rcpp::sourceCpp('src/updatePPLCpp.cpp')
data(meuse.grid)
candi <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(1:nrow(candi), candi), ncol = 3)
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
x.min <- 40
y.min <- 40
cutoff <- sqrt((x.max * x.max) + (y.max * y.max))
points <- 100
lags <- 7
lags.base <- 2
criterion <- "distribution"
lags.type <- "exponential"
pairs <- FALSE
# random selection
set.seed(2001)
res <- countPPL(points = points, candi = candi, lags = lags, 
                lags.type = lags.type, lags.base = lags.base, cutoff = cutoff,
                pairs = pairs)
set.seed(2001)
objPPL(points = points, candi = candi, lags = lags, lags.type = lags.type, 
       pairs = pairs, lags.base = lags.base, cutoff = cutoff, 
       criterion = criterion)
sum(points - res$points)
# vector of indexes
points <- 1:100
res <- countPPL(points = points, candi = candi, lags = lags, 
                lags.type = lags.type, lags.base = lags.base, cutoff = cutoff,
                pairs = pairs)
objPPL(points = points, candi = candi, lags = lags, lags.type = lags.type, 
       pairs = pairs, lags.base = lags.base, cutoff = cutoff, 
       criterion = criterion)
sum(length(points) - res$points)
#
# 3) Unit test #################################################################
rm(list = ls())
gc()
source('R/spSANNtools.R')
source('R/spJitter.R')
source('R/optimPPL.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
Rcpp::sourceCpp('src/updatePPLCpp.cpp')
data(meuse.grid)
candi <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(1:nrow(candi), candi), ncol = 3)
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
x.min <- 40
y.min <- 40
cutoff <- sqrt((x.max * x.max) + (y.max * y.max))
points <- 100
lags <- 1
lags.base <- 2
criterion <- "distribution"
lags.type <- "exponential"
pairs <- FALSE
set.seed(2001)
countPPL(points = points, candi = candi, lags = lags, lags.type = lags.type,
         lags.base = lags.base, cutoff = cutoff, pairs = pairs)
pairs <- TRUE
set.seed(2001)
countPPL(points = points, candi = candi, lags = lags, lags.type = lags.type,
         lags.base = lags.base, cutoff = cutoff, pairs = pairs)
