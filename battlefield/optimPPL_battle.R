# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(sp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)

# 0) DEFAULT EXAMPLE ###########################################################
require(sp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
schedule <- scheduleSPSANN(
  chains = 1,
  initial.acceptance = c(0.9, 0.99),
  initial.temperature = 6,
  x.max = 1540, y.max = 2060, x.min = 0,
  y.min = 0, cellsize = 40)
set.seed(2001)
res <- optimPPL(points = 10, candi = candi, schedule = schedule)
objSPSANN(res) - objPPL(points = res, candi = candi)
countPPL(points = res, candi = candi)

# 1) Point pairs with many chains #############################################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid, package = "sp")
candi <- meuse.grid[, 1:2]
schedule <- scheduleSPSANN(initial.temperature = 500)
set.seed(2001)
res <- optimPPL(
  points = 30, candi = candi, pairs = TRUE, schedule = schedule, plotit = TRUE)
objSPSANN(res)
objPPL(points = res, pairs = TRUE, candi = candi)
countPPL(points = res, candi = candi, pairs = TRUE)
plot(res)

# 2) Points per lag - select sample points from candi #########################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]

# random selection
points <- 100
set.seed(2001)
res <- countPPL(points = points, candi = candi, cutoff = 1000)
set.seed(2001)
objPPL(points = points, candi = candi, cutoff = 1000)
sum(points - res$ppl) # 346
# vector of indexes
points <- 1:100
res <- countPPL(points = points, candi = candi, cutoff = 1000)
objPPL(points = points, candi = candi, cutoff = 1000)
sum(length(points) - res$ppl) # 266

# 3) Unit test ################################################################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
set.seed(2001)
countPPL(points = 100, candi = candi, lags = 1, cutoff = Inf)[3] - 100
set.seed(2001)
countPPL(points = 100, candi = candi, lags = 1, pairs = TRUE, cutoff = Inf)[3] - 100 * 99 / 2

# 4) ADD TEN POINTS TO EXISTING SPATIAL SAMPLE ################################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
schedule <- scheduleSPSANN(
  chains = 500, stopping = 100, initial.temperature = 200, x.max = 1540, y.max = 2060, x.min = 0, y.min = 0,
  cellsize = 40)
free <- 10
set.seed(1984)
fixed <- candi[sample(1:nrow(candi), 30), ]
set.seed(2001)
res <- optimPPL(
  points = list(free = free, fixed = fixed), candi = candi, schedule = schedule, plotit = TRUE)
objSPSANN(res) - objPPL(points = res, candi = candi)
countPPL(points = res, candi = candi)
plot(res)
