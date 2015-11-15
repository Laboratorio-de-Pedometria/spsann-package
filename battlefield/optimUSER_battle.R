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
schedule <- scheduleSPSANN(chains = 1, initial.temperature = 30)

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
resUSER <- optimUSER(points = 100, fun = objUSER, lags = lags, n_lags = 9,
                     n_pts = 100, candi = candi, schedule = schedule)
timeUSER <- Sys.time() - timeUSER

# Run the optimization using the respective function implemented in spsann
set.seed(2001)
timePPL <- Sys.time()
resPPL <- optimPPL(points = 100, candi = candi, lags = lags, 
                   schedule = schedule)
timePPL <- Sys.time() - timePPL

# Compare results
timeUSER
timePPL
lapply(list(resUSER@points, resPPL@points), countPPL, candi = candi, 
       lags = lags)
objSPSANN(resUSER) - objSPSANN(resPPL)
