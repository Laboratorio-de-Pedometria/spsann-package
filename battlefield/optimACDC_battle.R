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
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
set.seed(2001)
res <- optimACDC(points = 100, candi = candi, covars = meuse.grid[, 5],
                 use.coords = TRUE, x.max = x.max, x.min = 40, y.max = y.max, 
                 y.min = 40, boundary = boundary, iterations = 1000,
                 nadir = list(sim = 10, save.sim = FALSE),
                 utopia = list(user = list(DIST = 0, CORR = 0)))
tail(attr(res, "energy")$obj, 1) # 0.4087331
objACDC(points = res, candi = candi, covars = meuse.grid[, 5],
        use.coords = TRUE, nadir = list(sim = 10, save.sim = FALSE), 
        utopia = list(user = list(DIST = 0, CORR = 0))) # 0.411831

# 1) FACTOR COVARIATES USING THE COORDINATES, WITH USER-DEFINED NADIR ##########
rm(list = ls())
gc()
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
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
nadir <- list(user = list(DIST = 10, CORR = 1))
utopia <- list(user = list(DIST = 0, CORR = 0))
set.seed(2001)
tmp <- optimACDC(points = 100, candi = candi, covars = meuse.grid[, 6:7], 
                 use.coords = TRUE, x.max = x.max, x.min = 40, 
                 y.max = y.max, y.min = 40, boundary = boundary, 
                 iterations = 100, nadir = nadir, utopia = utopia)
tail(attr(tmp, "energy")$obj, 1) # 1.724632
objACDC(points = tmp, candi = candi, covars = meuse.grid[, 6:7], 
        use.coords = TRUE, nadir = nadir, utopia = utopia)

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
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
nadir <- list(sim = 10, save.sim = TRUE)
utopia <- list(user = list(CORR = 0, DIST = 0))
set.seed(2001)
tmp <- optimACDC(points = 10, candi = candi, covars = meuse.grid[, 6:7],
                 use.coords = TRUE, x.max = x.max, x.min = 40, y.max = y.max, 
                 y.min = 40, boundary = boundary, iterations = 100, 
                 nadir = nadir, utopia = utopia)
tail(attr(tmp, "energy")$obj, 1) # 2.084277
objACDC(points = tmp, candi = candi, covars = meuse.grid[, 6:7], 
        use.coords = TRUE, nadir = nadir, utopia = utopia)

# 4) CATEGORICAL COVARIATES WITH MANY COVARIATES AND MANY POINTS ###############
# The function table() in the functions cramer() and chisqTest() is the one
# taking most of the time to run. Can we implement it in C++?
rm(list = ls())
gc()
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
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
set.seed(2001)
tmp <- optimACDC(points = 500, candi = candi, 
                 covars = meuse.grid[, rep(c(6, 7), 10)], 
                 use.coords = TRUE, x.max = x.max, x.min = 40, y.max = y.max, 
                 y.min = 40, boundary = boundary, iterations = 100, 
                 nadir = list(sim = 10, save.sim = TRUE), 
                 utopia = list(user = list(CORR = 0, DIST = 0)))
tail(attr(tmp, "energy")$obj, 1) # 4.75408
objACDC(points = tmp, candi = candi, covars = meuse.grid[, rep(c(6, 7), 10)],
        use.coords = TRUE, nadir = list(sim = 10, save.sim = TRUE), 
        utopia = list(user = list(CORR = 0, DIST = 0)))
