# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
require(Hmisc)
source('R/optimACDC.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
# 0) DEFAULT EXAMPLE ###########################################################
data(meuse.grid)
candi <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(1:dim(candi)[1], candi), ncol = 3)
str(meuse.grid)
covars <- meuse.grid[, 5]
x.max <- diff(bbox(boundary)[1, ])
y.min <- x.min <- 40
y.max <- diff(bbox(boundary)[2, ])
nadir <- list(sim = 10, save.sim = TRUE, user = NULL, abs = NULL)
weights <- list(strata = 0.5, correl = 0.5)
set.seed(2001)
res <- optimACDC(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, covars.type = "numeric", weights = weights,
                 x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min,
                 boundary = boundary, nadir = nadir, iterations = 100)
str(res)
#
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
nadir <- list(100)
set.seed(2001)
tmp <- optimACDC(points = points, candidates = candidates, covars = covars, continuous = continuous, use.coords = use.coords, x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, boundary = boundary, iterations = iterations, weights = weights, acceptance = acceptance, nadir = nadir)
# 2) CATEGORICAL COVARIATES ####################################################
covars <- meuse.grid[, c(6, 7)]
points <- 100
x.max <- diff(bbox(boundary)[1, ])
y.min <- x.min <- 40
y.max <- diff(bbox(boundary)[2, ])
iterations <- 500
acceptance <- list(initial = 0.99, cooling = iterations / 10)
weights <- list(strata = 0.5, correl = 0.5)
continuous <- FALSE
use.coords <- FALSE
sim.nadir <- 100
set.seed(2001)
tmp <- optimACDC(points = points, candidates = candidates, covars = covars, continuous = continuous, use.coords = use.coords, x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, boundary = boundary, iterations = iterations, weights = weights, acceptance = acceptance, sim.nadir = sim.nadir)
# 3) CATEGORICAL COVARIATES USING THE COORDINATES ##############################
# The following error appeared when the number of points is small (n = 5, 
# seed = 2001):
# Error in chisq.test(x[, i], x[, j], correct = FALSE) : 
#  'x' and 'y' must have at least 2 levels
# This error occurs because all points lie in the same class for one or more
# covariates.
covars <- meuse.grid[, c(6, 7)]
points <- 100
x.max <- diff(bbox(boundary)[1, ])
y.min <- x.min <- 40
y.max <- diff(bbox(boundary)[2, ])
iterations <- 500
acceptance <- list(initial = 0.99, cooling = iterations / 10)
weights <- list(strata = 0.5, correl = 0.5)
continuous <- FALSE
use.coords <- TRUE
stopping <- list(max.count = iterations / 10)
sim.nadir <- 10
strata.type <- "equal.area"
plotit <- TRUE
progress <- TRUE
verbose <- TRUE
set.seed(2001)
tmp <- optimACDC(points = points, candidates = candidates, covars = covars, continuous = continuous, use.coords = use.coords, x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, boundary = boundary, iterations = iterations, weights = weights, acceptance = acceptance, sim.nadir = sim.nadir, stopping = stopping, plotit = plotit, progress = progress, verbose = verbose, strata.type = strata.type)
# 4) CATEGORICAL COVARIATES WITH MANY COVARIATES AND MANY POINTS ###############
# The function table() in the functions cramer() and chisqTest() is the one
# taking most of the time to run. Can we implement it in C++?
covars <- meuse.grid[, rep(c(6, 7), 10)]
points <- 500
x.max <- diff(bbox(boundary)[1, ])
y.min <- x.min <- 40
y.max <- diff(bbox(boundary)[2, ])
iterations <- 2
acceptance <- list(initial = 0.99, cooling = iterations / 10)
weights <- list(strata = 0.5, correl = 0.5)
continuous <- FALSE
use.coords <- FALSE
stopping <- list(max.count = iterations / 10)
sim.nadir <- 1
strata.type <- "equal.area"
plotit <- TRUE
progress <- TRUE
verbose <- FALSE
set.seed(2001)
Rprof(a <- tempfile())
tmp <- optimACDC(points = points, candidates = candidates, covars = covars, continuous = continuous, use.coords = use.coords, x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, boundary = boundary, iterations = iterations, weights = weights, acceptance = acceptance, sim.nadir = sim.nadir, stopping = stopping, plotit = plotit, progress = progress, verbose = verbose, strata.type = strata.type)
Rprof()
summaryRprof(a)
unlink(a)
