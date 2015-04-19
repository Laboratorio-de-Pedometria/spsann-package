# Initial settings
rm(list = ls())
gc()
source('R/spSANNtools.R')
source('R/spJitter.R')
source('R/optimUSER.R')
source('R/optimPPL.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
Rcpp::sourceCpp('src/updatePPLCpp.cpp')
# 0) DEFAULT EXAMPLE ###########################################################
require(pedometrics)
require(sp)
require(rgeos)
require(SpatialTools)
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

# Define the objective function - number of points per lag distance class
objUSER <-
  function (points, lags, n_lags, n_pts) {
  dm <- SpatialTools::dist1(points[, 2:3])
  ppl <- vector()
  for (i in 1:n_lags) {
    n <- which(dm > lags[i] & dm <= lags[i + 1], arr.ind = TRUE)
    ppl[i] <- length(unique(c(n)))
  }
  distri <- rep(n_pts, n_lags)
  res <- sum(distri - ppl)
}
lags <- seq(1, 1000, length.out = 10)

# Run the optimization using the user-defined objective function
set.seed(2001)
timeUSER <- Sys.time()
resUSER <- optimUSER(points = 100, fun = objUSER, lags = lags, 
                     n_lags = 9, n_pts = 100,
                     candi = candi, x.max = x.max, x.min = 40, y.max = y.max,
                     y.min = 40, boundary = boundary, iterations = 1000)
timeUSER <- Sys.time() - timeUSER

# Run the optimization using the respective function implemented in spsann
set.seed(2001)
timePPL <- Sys.time()
resPPL <- optimPPL(points = 100, candi = candi, lags = lags,  
                   criterion = "distribution", x.max = x.max, x.min = 40, 
                   y.max = y.max, y.min = 40, boundary = boundary,
                   iterations = 1000)
timePPL <- Sys.time() - timePPL

# Compare results
timeUSER
timePPL
lapply(list(resUSER, resPPL), countPPL, lags = lags, pairs = FALSE)
x <- attr(resUSER, "energy.state")
y <- attr(resPPL, "energy.state")
sapply(list(x, y), tail, 1)
plot(x, y, asp = 1)
abline(0, 1, col = "red")
