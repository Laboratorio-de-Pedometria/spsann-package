# Initial settings
rm(list = ls())
gc()
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
# 0) DEFAULT EXAMPLE ###########################################################
require(sp)
require(SpatialTools)
data(meuse.grid)
candi <- meuse.grid[, 1:2]

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
                     n_pts = 100, candi = candi, iterations = 100,
                     plotit = FALSE, track = FALSE, verbose = FALSE)
timeUSER <- Sys.time() - timeUSER

# Run the optimization using the respective function implemented in spsann
set.seed(2001)
timePPL <- Sys.time()
resPPL <- optimPPL(points = 100, candi = candi, lags = lags, 
                   iterations = 100, plotit = FALSE, track = FALSE, 
                   verbose = FALSE)
timePPL <- Sys.time() - timePPL

# Compare results
timeUSER
timePPL
lapply(list(resUSER, resPPL), countPPL, lags = lags, pairs = FALSE)
x <- attr(resUSER, "energy.state") # 58
y <- attr(resPPL, "energy.state") # 58
sapply(list(x, y), tail, 1)
plot(x, y, asp = 1)
abline(0, 1, col = "red")
