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
cutoff <- sqrt((x.max * x.max) + (y.max * y.max)) / 2
set.seed(2001)
res <- optimPPL(points = 100, candi = candi, lags = 7, pairs = FALSE,
                lags.base = 2, criterion = "distribution", cutoff = cutoff,
                lags.type = "exponential", x.max = x.max, x.min = 40, 
                y.max = y.max, y.min = 40, boundary = boundary,
                iterations = 100, plotit = TRUE, verbose = TRUE)
countPPL(points = res, lags = 7, lags.type = "exponential", pairs = FALSE,
         lags.base = 2, cutoff = cutoff)
tail(attr(res, "energy.state"), 1) # 168
objPPL(points = res, lags = 7, lags.type = "exponential", pairs = FALSE,
       lags.base = 2, cutoff = cutoff, criterion = "distribution")

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
tail(attr(res, "energy.state"), 1) # 7608.857
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
candi <- matrix(cbind(1:nrow(candi), candi), ncol = 3)
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
cutoff <- sqrt((x.max * x.max) + (y.max * y.max)) / 2
set.seed(2001)
res <- optimPPL(points = 100, candi = candi, lags = 7, pairs = FALSE,
                lags.base = 2, criterion = "distribution", cutoff = cutoff,
                lags.type = "exponential", x.max = x.max, x.min = 40, 
                y.max = y.max, y.min = 40, boundary = boundary,
                iterations = 100, plotit = TRUE, verbose = TRUE, greedy = TRUE)
countPPL(points = res, lags = 7, lags.type = "exponential", pairs = FALSE,
         lags.base = 2, cutoff = cutoff)
tail(attr(res, "energy.state"), 1) # 162
objPPL(points = res, lags = 7, lags.type = "exponential", pairs = FALSE,
       lags.base = 2, cutoff = cutoff, criterion = "distribution")

# 5) OLD INFINITE LOOP #########################################################
rm(list = ls())
gc()
source('R/spSANNtools.R')
source('R/spJitter.R')
source('R/optimPPL.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
Rcpp::sourceCpp('src/updatePPLCpp.cpp')
nx <- ny <- 100
s1 <- 1:nx
s1 <- s1 - 0.5
s2 <- s1
candi <- expand.grid(s1 = s1, s2 = s2)
coordinates(candi) <- ~ s1 + s2
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(1:nrow(candi), candi), ncol = 3)
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
cutoff <- sqrt((x.max * x.max) + (y.max * y.max))# / 2
lags <- 7
points <- 49
set.seed(2001)
sample_b <- optimPPL(points = points, candi = candi, cutoff = cutoff, 
                     lags = lags, x.max = x.max, x.min = 2, y.max = y.max, 
                     y.min = 2, boundary = boundary, iterations = 7000, 
                     plotit = FALSE, pairs = TRUE, verbose = FALSE, 
                     greedy = TRUE)
countPPL(sample_b, lags = lags, cutoff = cutoff, pairs = TRUE)
objPPL(sample_b, lags = lags, cutoff = cutoff, pairs = TRUE)
tail(attr(sample_b, "energy.state"), 1) # 1176
