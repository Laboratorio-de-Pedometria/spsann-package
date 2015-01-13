# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
require(Hmisc)
source('R/spSANNtools.R')
source('R/spJitter.R')
source('R/optimACDC.R')
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
# PPL
lags.type <- "exponential"
ppl.criterion <- "distribution"
lags <- 7
lags.base <- 2
pre.distri <- NULL
cutoff <- sqrt((x.max * x.max) + (y.max * y.max))
# ACDC
covars <- meuse.grid[, 1]
covar.type <- "numeric"
use.coords <- TRUE
strata.type <- "equal.area"
weights.ACDC <- list(strata = 0.5, correl = 0.5)
# MOOP
weights.PAN <- list(PPL = 1/3, ACDC = 1/3, MSSD = 1/3)
nadir <- list(sim = 1000, user = NULL, abs = NULL)
