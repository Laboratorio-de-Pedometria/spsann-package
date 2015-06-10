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
require(SpatialTools)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
set.seed(2001)
res <- optimPPL(points = 100, candi = candi, iterations = 100)
#countPPL(points = res, lags = 7, lags.type = "exponential", pairs = FALSE,
#         lags.base = 2, cutoff = cutoff)
tail(attr(res, "energy.state"), 1) # 160
#objPPL(points = res, lags = 7, lags.type = "exponential", pairs = FALSE,
#       lags.base = 2, cutoff = cutoff, criterion = "distribution")
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
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
cutoff <- sqrt((x.max * x.max) + (y.max * y.max))
pairs <- TRUE
set.seed(2001)
res <- optimPPL(points = 100, candi = candi, pairs = pairs, cutoff = cutoff,
                x.max = x.max, x.min = 40, y.max = y.max, y.min = 40, 
                boundary = boundary, iterations = 100)
countPPL(points = res, cutoff = cutoff, pairs = pairs)
tail(attr(res, "energy.state"), 1) # 7608.857
objPPL(points = res, pairs = pairs, cutoff = cutoff)

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
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
cutoff <- sqrt((x.max * x.max) + (y.max * y.max))
points <- 100
# random selection
set.seed(2001)
res <- countPPL(points = points, candi = candi, cutoff = cutoff)
set.seed(2001)
objPPL(points = points, candi = candi, cutoff = cutoff)
sum(points - res$points) # 136
# vector of indexes
points <- 1:100
res <- countPPL(points = points, candi = candi, cutoff = cutoff)
objPPL(points = points, candi = candi, cutoff = cutoff)
sum(length(points) - res$points) # 300

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
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
cutoff <- sqrt((x.max * x.max) + (y.max * y.max))
set.seed(2001)
countPPL(points = 100, candi = candi, lags = 1, cutoff = cutoff)
set.seed(2001)
countPPL(points = 100, candi = candi, lags = 1, cutoff = cutoff, pairs = TRUE)

# 4) GREEDY ALGORITH ###########################################################
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
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
cutoff <- sqrt((x.max * x.max) + (y.max * y.max)) / 2
set.seed(2001)
res <- optimPPL(points = 100, candi = candi, cutoff = cutoff, x.max = x.max,
                x.min = 40, y.max = y.max, y.min = 40, boundary = boundary,
                iterations = 100, greedy = TRUE)
countPPL(points = res, cutoff = cutoff)
tail(attr(res, "energy.state"), 1) # 162
objPPL(points = res, cutoff = cutoff)

# 5) OLD INFINITE LOOP #########################################################
rm(list = ls())
gc()
source('R/spSANNtools.R')
source('R/spJitter.R')
source('R/optimPPL.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
Rcpp::sourceCpp('src/updatePPLCpp.cpp')
nx <- ny <- 100
x <- 1:nx
x <- x - 0.5
y <- x
candi <- expand.grid(x = x, y = y)
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
cutoff <- sqrt((x.max * x.max) + (y.max * y.max))# / 2
points <- 49
set.seed(2001)
sample_b <- optimPPL(points = points, candi = candi, cutoff = cutoff, 
                     x.max = x.max, x.min = 2, y.max = y.max, 
                     y.min = 2, boundary = boundary, iterations = 7000, 
                     plotit = FALSE, pairs = TRUE, verbose = FALSE, 
                     greedy = TRUE)
countPPL(sample_b, cutoff = cutoff, pairs = TRUE)
tail(attr(sample_b, "energy.state"), 1) # 1176
objPPL(sample_b, cutoff = cutoff, pairs = TRUE)

