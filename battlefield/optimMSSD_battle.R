# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(SpatialTools)
require(sp)
require(rgeos)
source('R/spSANNtools.R')
source('R/spJitter.R')
source('R/optimMSSD.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
Rcpp::sourceCpp('src/calcMSSDCpp.cpp')
Rcpp::sourceCpp('src/updateMSSDCpp.cpp')
# PREPARE DATA #################################################################
data(meuse.grid)
candidates <- meuse.grid[, 1:2]
coordinates(candidates) <- ~ x + y
gridded(candidates) <- TRUE
boundary <- as(candidates, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candidates <- coordinates(candidates)
candidates <- matrix(cbind(c(1:dim(candidates)[1]), candidates), ncol = 3)
points <- 100
x.max <- diff(bbox(boundary)[1, ])
y.min <- x.min <- 40
y.max <- diff(bbox(boundary)[2, ])
iterations <- 10000
acceptance <- list(initial = 0.99, cooling = iterations / 10)
weights <- list(strata = 0.5, correl = 0.5)
stopping <- list(max.count = iterations / 10)
plotit <- TRUE
progress <- TRUE
verbose <- TRUE
# 1) FIRST TEST ################################################################
x11()
set.seed(2001)
tmp <- optimMSSD(points = points, candidates = candidates, x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, iterations = iterations, acceptance = acceptance, stopping = stopping, plotit = plotit, boundary = boundary, progress = progress, verbose = verbose)
