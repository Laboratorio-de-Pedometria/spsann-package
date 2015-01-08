# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(sp)
require(rgeos)
source('~/PROJECTS/r-packages/spsann/R/optimACDC.R')
source('~/PROJECTS/r-packages/spsann/R/spSANNtools.R')
source('~/PROJECTS/r-packages/pedometrics/cooking/utils.R')
source('~/PROJECTS/r-packages/pedometrics/pkg/pedometrics/R/cramer.R')
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
# 1) CONTINUOUS COVARIATES USIGN THE COORDINATES ###############################
covars <- meuse.grid[, 5]
points <- 100
x.max <- diff(bbox(boundary)[1, ])
y.min <- x.min <- 40
y.max <- diff(bbox(boundary)[2, ])
iterations <- 100
acceptance <- list(initial = 0.99, cooling = iterations / 10)
weights <- list(strata = 0.5, correl = 0.5)
continuous <- TRUE
use.coords <- TRUE
sim.nadir <- 100
set.seed(2001)
tmp <- optimACDC(points = points, candidates = candidates, covars = covars, 
                 continuous = continuous, use.coords = use.coords, 
                 x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, 
                 boundary = boundary, iterations = iterations, 
                 weights = weights, acceptance = acceptance, 
                 sim.nadir = sim.nadir)

# 2) CATEGORICAL COVARIATES ####################################################
covars <- meuse.grid[, c(6, 7)]
points <- 100
x.max <- diff(bbox(boundary)[1, ])
y.min <- x.min <- 40
y.max <- diff(bbox(boundary)[2, ])
iterations <- 1000
acceptance <- list(initial = 0.99, cooling = iterations / 10)
weights <- list(strata = 0.5, correl = 0.5)
continuous <- FALSE
use.coords <- FALSE
sim.nadir <- 100
set.seed(2001)
tmp <- optimACDC(points = points, candidates = candidates, covars = covars, 
                 continuous = continuous, use.coords = use.coords, 
                 x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, 
                 boundary = boundary, iterations = iterations, 
                 weights = weights, acceptance = acceptance, 
                 sim.nadir = sim.nadir)
#
# 2) CATEGORICAL COVARIATES USING THE COORDINATES ##############################
# Error in breaks[, i] : incorrect number of dimensions
#
cont2cat <-
  function (covars, breaks) {
    for (i in 1:ncol(covars)) {
      covars[, i] <- cut2(covars[, i], breaks[[i]])
    }
    return (covars)
  }
covars <- meuse.grid[, c(6, 7)]
points <- 5
x.max <- diff(bbox(boundary)[1, ])
y.min <- x.min <- 40
y.max <- diff(bbox(boundary)[2, ])
iterations <- 1000
acceptance <- list(initial = 0.99, cooling = iterations / 10)
weights <- list(strata = 0.5, correl = 0.5)
continuous <- FALSE
use.coords <- TRUE
stopping <- list(max.count = iterations / 10)
sim.nadir <- 1000
strata.type <- "equal.area"
plotit <- TRUE
progress <- TRUE
verbose <- TRUE
set.seed(2001)
tmp <- optimACDC(points = points, candidates = candidates, covars = covars, 
                 continuous = continuous, use.coords = use.coords, 
                 x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, 
                 boundary = boundary, iterations = iterations, 
                 weights = weights, acceptance = acceptance, 
                 sim.nadir = sim.nadir, stopping = stopping, plotit = plotit,
                 progress = progress, verbose = verbose, 
                 strata.type = strata.type)
