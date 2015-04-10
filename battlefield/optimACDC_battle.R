# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
source('R/optimACDC.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
# 0) DEFAULT EXAMPLE ###########################################################
require(pedometrics)
require(sp)
require(rgeos)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(1:dim(candi)[1], candi), ncol = 3)
covars <- meuse.grid[, 5]
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
nadir <- list(sim = 10, save.sim = TRUE, user = NULL, abs = NULL)
utopia <- list(user = list(correl = 0, strata = 0), abs = NULL)
scale <- list(type = "upper-lower", max = 100)
weights <- list(strata = 0.5, correl = 0.5)
set.seed(2001)
res <- optimACDC(points = 100, candi = candi, covars = covars,
                 use.coords = TRUE, weights = weights, x.max = x.max, 
                 x.min = 40, y.max = y.max, y.min = 40,
                 boundary = boundary, nadir = nadir, iterations = 1000,
                 utopia = utopia, scale = scale)
tail(attr(res, "energy"), 1) # 55.59217
objACDC(points = res, candi = candi, covars = covars, use.coords = TRUE, 
        weights = weights, nadir = nadir, utopia = utopia, scale = scale)

# 1) FACTOR COVARIATES USING THE COORDINATES ###################################
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
source('R/optimACDC.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
data(meuse.grid)
candi <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(1:dim(candi)[1], candi), ncol = 3)
covars <- meuse.grid[, 6:7]
covars.type <- "factor"
points <- 100
x.max <- diff(bbox(boundary)[1, ])
y.min <- x.min <- 40
y.max <- diff(bbox(boundary)[2, ])
iterations <- 100
acceptance <- list(initial = 0.99, cooling = iterations / 10)
weights <- list(strata = 0.5, correl = 0.5)
use.coords <- TRUE
nadir <- list(sim = 10, save.sim = TRUE, user = NULL, abs = NULL)
utopia <- list(user = list(correl = 0, strata = 0), abs = NULL)
scale <- list(type = "upper-lower", max = 100)
set.seed(2001)
tmp <- optimACDC(points = points, candi = candi, covars = covars, 
                 covars.type = covars.type, use.coords = use.coords, 
                 x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min,
                 boundary = boundary, iterations = iterations, 
                 weights = weights, acceptance = acceptance, nadir = nadir,
                 utopia = utopia, scale = scale)
tail(attr(tmp, "energy"), 1) # 60.25659
objACDC(points = tmp, candi = candi, covars = covars, use.coords = use.coords, 
        covars.type = covars.type, weights = weights, nadir = nadir,
        utopia = utopia, scale = scale) # 60.25659
# 3) FACTOR COVARIATES USING THE COORDINATES WITH A FEW POINTS #################
# The following error appeared in an old version (before correcting for the 
# number of strata) when the number of points is small (n = 5, seed = 2001):
# Error in chisq.test(x[, i], x[, j], correct = FALSE) : 
#  'x' and 'y' must have at least 2 levels
# This error seems to occur because all points lie in the same class for one 
# or more covariates.
# ERROR: The following error appears when the number of points is small (n < 10,
# seed = 2001):
# Error in if (new_energy <= old_energy) { : 
# missing value where TRUE/FALSE needed
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
source('R/optimACDC.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
data(meuse.grid)
candi <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
candi <- matrix(cbind(1:dim(candi)[1], candi), ncol = 3)
covars <- meuse.grid[, 6:7]
covars.type <- "factor"
use.coords <- TRUE
x.max <- diff(bbox(boundary)[1, ])
y.min <- x.min <- 40
y.max <- diff(bbox(boundary)[2, ])
iterations <- 100
acceptance <- list(initial = 0.99, cooling = iterations / 10)
weights <- list(strata = 0.5, correl = 0.5)
nadir <- list(sim = 10, save.sim = TRUE, user = NULL, abs = NULL)
utopia <- list(user = list(correl = 0, strata = 0), abs = NULL)
scale <- list(type = "upper-lower", max = 100)
points <- 10
set.seed(2001)
tmp <- optimACDC(points = points, candi = candi, covars = covars,
                 covars.type = covars.type, use.coords = use.coords,
                 x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min,
                 boundary = boundary, iterations = iterations,
                 weights = weights, acceptance = acceptance, nadir = nadir,
                 utopia = utopia, scale = scale)
objACDC(points = points, candi = candi, covars = covars,
        covars.type = covars.type, use.coords = use.coords, weights = weights,
        nadir = nadir, utopia = utopia, scale = scale)
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
