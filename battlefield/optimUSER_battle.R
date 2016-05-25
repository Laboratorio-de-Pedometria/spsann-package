# Initial settings
rm(list = ls())
gc()
require(pedometrics)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)

# 0) DEFAULT EXAMPLE ###########################################################
require(sp)
require(SpatialTools)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
schedule <- scheduleSPSANN(chains = 1, initial.temperature = 30,
                           x.max = 1540, y.max = 2060, x.min = 0, 
                           y.min = 0, cellsize = 40)

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
resUSER <- optimUSER(points = 10, fun = objUSER, lags = lags, n_lags = 9,
                     n_pts = 10, candi = candi, schedule = schedule)
timeUSER <- Sys.time() - timeUSER

# Run the optimization using the respective function implemented in spsann
set.seed(2001)
timePPL <- Sys.time()
resPPL <- optimPPL(points = 10, candi = candi, lags = lags, 
                   schedule = schedule)
timePPL <- Sys.time() - timePPL

# Compare results
timeUSER
timePPL
lapply(list(resUSER, resPPL), countPPL, candi = candi, lags = lags)
objSPSANN(resUSER) - objSPSANN(resPPL)

# 1) Error found by Alexandre Wadoux on 25 May 2016 ###########################################################
# Message:
# I am running my final processing with the new optimUSER function in spsann. But I don’t understand why the
# process stops after each chain. For example if I put chains=2 then the optimisation process stops a 50% and 
# returns the object. Is it normal ? Should the scheduleSPSANN object be done within a loop ?
# # I don’t find this problem in the package documentation this is why I am wondering if it is normal.

rm(list = ls())
gc()
data(meuse.grid, package = "sp")
candi <- meuse.grid[, 1:2]
schedule <- spsann::scheduleSPSANN(
  chains = 2, initial.temperature = 30, x.max = 1540, y.max = 2060, x.min = 0, y.min = 0, cellsize = 40)

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
resUSER <- spsann::optimUSER(
  points = 100, fun = objUSER, lags = lags, n_lags = 9, n_pts = 10, candi = candi, schedule = schedule, 
  track = TRUE, plotit = TRUE)

# Check results
# If the error occurs, then the number of recorded energy values shoudl be 0.5 * chains = points
nrow(resUSER$objective$energy)
