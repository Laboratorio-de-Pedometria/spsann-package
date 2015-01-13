# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
require(Hmisc)
source('~/PROJECTS/r-packages/spsann/R/spSANNtools.R')
source('~/PROJECTS/r-packages/spsann/R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
Rcpp::sourceCpp('src/updatePPLCpp.cpp')
# PREPARE DATA #################################################################
data(meuse.grid)
candidates <- meuse.grid[, 1:2]
coordinates(candidates) <- ~ x + y
gridded(candidates) <- TRUE
boundary <- as(candidates, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candidates <- coordinates(candidates)
candidates <- matrix(cbind(c(1:dim(candidates)[1]), candidates), ncol = 3)
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
tmp <- optimPPL(points = points, candidates = candidates, lags = lags, lags.type = lags.type, lags.base = lags.base, cutoff = cutoff, criterion = criterion, pre.distri = pre.distri, x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, boundary = boundary, iterations = iterations, acceptance = acceptance, stopping = stopping)
pointsPerLag(points = tmp, lags = lags, lags.type = lags.type, cutoff = cutoff)
objPoints(points = tmp, lags = lags, lags.type = lags.type, lags.base = lags.base, cutoff = cutoff, criterion = criterion, pre.distri = pre.distri)
