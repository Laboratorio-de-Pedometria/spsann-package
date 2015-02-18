# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
source('R/optimDIST.R')
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
y.min <- 40
x.min <- 40
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, covars.type = "numeric", x.max = x.max, 
                 x.min = x.min, y.max = y.max, y.min = y.min,
                 boundary = boundary, iterations = 500)
tail(attr(res, "energy"))
objDIST(points = res, candi = candi, covars = covars, covars.type = "numeric",
        use.coords = TRUE)
# 1) FACTOR COVARIATES #########################################################
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
source('R/optimCORR.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
data(meuse.grid)
candi              <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi)     <- TRUE
boundary           <- as(candi, "SpatialPolygons")
boundary           <- gUnionCascaded(boundary)
candi              <- coordinates(candi)
candi              <- matrix(cbind(1:dim(candi)[1], candi), ncol = 3)
covars             <- meuse.grid[, 6:7]
x.max              <- diff(bbox(boundary)[1, ])
y.max              <- diff(bbox(boundary)[2, ])
y.min              <- 40
x.min              <- 40
set.seed(2001)
res <- optimCORR(points = 100, candi = candi, covars = covars, 
                 use.coords = FALSE, covars.type = "factor", x.max = x.max, 
                 x.min = x.min, y.max = y.max, y.min = y.min,
                 boundary = boundary, iterations = 500)
tail(attr(res, "energy"))
objCORR(points = res, candi = candi, covars = covars, covars.type = "factor",
        use.coords = FALSE)
# 2) FACTOR COVARIATES WITH THE COORDINATES ####################################
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
source('R/optimCORR.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
data(meuse.grid)
candi              <- meuse.grid[, 1:2]
coordinates(candi) <- ~ x + y
gridded(candi)     <- TRUE
boundary           <- as(candi, "SpatialPolygons")
boundary           <- gUnionCascaded(boundary)
candi              <- coordinates(candi)
candi              <- matrix(cbind(1:dim(candi)[1], candi), ncol = 3)
covars             <- meuse.grid[, 6:7]
x.max              <- diff(bbox(boundary)[1, ])
y.max              <- diff(bbox(boundary)[2, ])
y.min              <- 40
x.min              <- 40
set.seed(2001)
res <- optimCORR(points = 100, candi = candi, covars = covars, 
                 strata.type = "area", use.coords = TRUE,
                 covars.type = "factor", x.max = x.max, x.min = x.min, 
                 y.max = y.max, y.min = y.min, 
                 boundary = boundary, iterations = 500)
tail(attr(res, "energy"))
objCORR(points = res, candi = candi, covars = covars, covars.type = "factor",
        use.coords = TRUE)
