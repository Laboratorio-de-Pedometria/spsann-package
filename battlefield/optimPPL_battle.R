# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(SpatialTools)
require(rgeos)
source('R/spSANNtools.R')
source('R/spJitter.R')
source('R/optimPPL.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
Rcpp::sourceCpp('src/updatePPLCpp.cpp')
# 0) DEFAULT EXAMPLE ###########################################################
# Check this! The optimization is returning an energy value different from that
# returned by the function designed to calculate the energy
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
cutoff <- sqrt((x.max * x.max) + (y.max * y.max))
points <- 100
set.seed(2002)
res <- optimPPL(points = points, candi = candi, lags = 7, lags.base = 2,
                criterion = "distribution", lags.type = "exponential",
                cutoff = cutoff, x.max = x.max, x.min = 40, y.max = y.max, 
                y.min = 40, boundary = boundary, iterations = 1000, 
                plotit = TRUE)
str(res)
pointsPerLag(points = res, lags = 7, lags.type = "exponential", lags.base = 2, 
             cutoff = cutoff)
objPoints(points = res, lags = 7, lags.type = "exponential", lags.base = 2,
          cutoff = cutoff, criterion = "distribution")
tail(attr(res, "energy.state"), 1)
# PREPARE DATA #################################################################
data(meuse.grid)
candi <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(1:nrow(candi), candi), ncol = 3)
#
# 1) DISTRIBUTION = EXPONENTIAL; CRITERION = DISTRIBUTION; POINTS = 100 ########
lags.type <- "exponential"
criterion <- "distribution"
points <- 100
lags <- 7
lags.base <- 2
pre.distri <- NULL
x.max <- diff(bbox(boundary)[1, ])
y.min <- x.min <- 40
y.max <- diff(bbox(boundary)[2, ])
cutoff <- sqrt((x.max * x.max) + (y.max * y.max))
iterations <- 10000
acceptance <- list(initial = 0.99, cooling = iterations / 10)
weights <- list(strata = 0.5, correl = 0.5)
stopping <- list(max.count = iterations / 10)
plotit <- TRUE
progress <- TRUE
verbose <- TRUE
X11()
set.seed(2001)
tmp <- optimPPL(points = points, candi = candi, lags = lags, lags.type = lags.type, lags.base = lags.base, cutoff = cutoff, criterion = criterion, pre.distri = pre.distri, x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, boundary = boundary, iterations = iterations, acceptance = acceptance, stopping = stopping)
pointsPerLag(points = tmp, lags = lags, lags.type = lags.type, cutoff = cutoff)
objPoints(points = tmp, lags = lags, lags.type = lags.type, lags.base = lags.base, cutoff = cutoff, criterion = criterion, pre.distri = pre.distri)
