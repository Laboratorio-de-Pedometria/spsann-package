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
nadir <- list(sim = 10, seeds = 1:10)
utopia <- list(user = list(DIST = 0, CORR = 0))
covars <- meuse.grid[, 5]
set.seed(2001)
res <- optimACDC(points = 100, candi = candi, covars = covars, y.max = y.max,
                 use.coords = TRUE, x.max = x.max, x.min = 40, y.min = 40, 
                 boundary = boundary, iterations = 100, nadir = nadir, 
                 utopia = utopia)
tail(attr(res, "energy")$obj, 1) # 0.5251647
objACDC(points = res, candi = candi, covars = covars, use.coords = TRUE, 
        nadir = nadir, utopia = utopia)
# MARGINAL DISTRIBUTION
par(mfrow = c(3, 3))
# Covariates
i <- sample(1:nrow(candi), 100)
hist(candi[, 1], breaks = 10)
hist(candi[, 2], breaks = 10)
hist(covars, breaks = 10)
# Optimized sample
hist(candi[res[, 1], 1], breaks = 10)
hist(candi[res[, 1], 2], breaks = 10)
hist(covars[res[, 1]], breaks = 10)
# Random sample
hist(candi[i, 1], breaks = 10)
hist(candi[i, 2], breaks = 10)
hist(covars[i], breaks = 10)

# LINEAR CORRELATION
# Covariates
cor(cbind(candi[, 1], candi[, 2], covars))
# Optimized sample
cor(cbind(candi[res[, 1], 1], candi[res[, 1], 2], covars[res[, 1]]))
# Random sample
cor(cbind(candi[i, 1], candi[i, 2], covars[i]))

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
tail(attr(tmp, "energy")$obj, 1) # 1.736775
objACDC(points = tmp, candi = candi, covars = meuse.grid[, 6:7], 
        use.coords = TRUE, nadir = nadir, utopia = utopia)

# 3) FACTOR COVARIATES USING THE COORDINATES WITH A FEW POINTS #################
# Tue 9 Jun: objACDC() does not return the same criterion value if 
#            'iterations = 100'
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
covars <- meuse.grid[, 6:7]
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
nadir <- list(sim = 10, seeds = 1:10)
utopia <- list(user = list(CORR = 0, DIST = 0))
set.seed(2001)
tmp <- optimACDC(points = 10, candi = candi, covars = covars,
                 use.coords = TRUE, x.max = x.max, x.min = 40, y.max = y.max, 
                 y.min = 40, boundary = boundary, iterations = 1000, 
                 nadir = nadir, utopia = utopia)
tail(attr(tmp, "energy")$obj, 1) # 0.7085472
objACDC(points = tmp, candi = candi, covars = covars, use.coords = TRUE,
        nadir = nadir, utopia = utopia)

# 4) CATEGORICAL COVARIATES WITH MANY COVARIATES AND MANY POINTS ###############
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
nadir <- list(sim = 10, seeds = 1:10)
set.seed(2001)
tmp <- optimACDC(points = 500, candi = candi, 
                 covars = meuse.grid[, rep(c(6, 7), 10)], 
                 use.coords = TRUE, x.max = x.max, x.min = 40, y.max = y.max, 
                 y.min = 40, boundary = boundary, iterations = 100, 
                 nadir = nadir, utopia = list(user = list(CORR = 0, DIST = 0)))
tail(attr(tmp, "energy")$obj, 1) # 0.5870467
objACDC(points = tmp, candi = candi, covars = meuse.grid[, rep(c(6, 7), 10)],
        use.coords = TRUE, nadir = nadir, 
        utopia = list(user = list(CORR = 0, DIST = 0)))

# 5) NUMERIC COVARIATES USING THE COORDINATES, WITH USER-DEFINED NADIR #########
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
covars = meuse.grid[, 5]
set.seed(2001)
tmp <- optimACDC(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, x.max = x.max, x.min = 40, 
                 y.max = y.max, y.min = 40, boundary = boundary, 
                 iterations = 100, nadir = nadir, utopia = utopia)
tail(attr(tmp, "energy")$obj, 1) # 0.1400659
objACDC(points = tmp, candi = candi, covars = covars, 
        use.coords = TRUE, nadir = nadir, utopia = utopia)
