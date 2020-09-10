# 0) DEFAULT EXAMPLE ##########################################################################################
rm(list = ls())
gc()
devtools::load_all()
# Example
data(meuse.grid, package = "sp")
candi <- meuse.grid[, 1:2]
schedule <- scheduleSPSANN(
  chains = 1,
  initial.acceptance = c(0.8, 0.99),
  initial.temperature = 9.5,
  x.max = 1540, y.max = 2060, x.min = 0,
  y.min = 0, cellsize = 40)
set.seed(2001)
res <- optimPPL(points = 10, candi = candi, schedule = schedule, track = TRUE)
objPPL(res)
countPPL(res)
plot(res)

# 1) Point pairs with many chains #############################################################################
rm(list = ls())
gc()
devtools::load_all()
data(meuse.grid, package = "sp")
candi <- meuse.grid[, 1:2]
set.seed(2001, sample.kind = "Round")
res <- optimPPL(
  points = 30, candi = candi, pairs = TRUE, plotit = TRUE,
  schedule = scheduleSPSANN(initial.temperature = 100, cellsize = 40))
objPPL(res)
countPPL(res)
plot(res)

# 2) Points per lag - select sample points from candi #########################################################
rm(list = ls())
gc()
devtools::load_all()
data(meuse.grid, package = "sp")
candi <- meuse.grid[, 1:2]
# random selection
points <- 100
set.seed(2001, sample.kind = "Round")
res <- countPPL(points = points, candi = candi, cutoff = 1000)
set.seed(2001, sample.kind = "Rounding")
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
devtools::load_all()
data(meuse.grid, package = "sp")
candi <- meuse.grid[, 1:2]
set.seed(2001, sample.kind = "Rounding")
countPPL(points = 100, candi = candi, lags = 1, cutoff = Inf)[3] - 100
set.seed(2001, sample.kind = "Rounding")
countPPL(points = 100, candi = candi, lags = 1, pairs = TRUE, cutoff = Inf)[3] - 100 * 99 / 2

# 4) ADD TEN POINTS TO EXISTING SPATIAL SAMPLE ################################################################
rm(list = ls())
gc()
devtools::load_all()
data(meuse.grid, package = "sp")
candi <- meuse.grid[, 1:2]
schedule <- scheduleSPSANN(
  initial.acceptance = c(0.9, 0.99), chains = 500, stopping = 100, initial.temperature = 120,
  x.max = 1540, y.max = 2060, x.min = 0, y.min = 0, cellsize = 40)
free <- 10
set.seed(1984, sample.kind = "Round")
fixed <- candi[sample(1:nrow(candi), 30), ]
set.seed(2001, sample.kind = "Round")
res <- optimPPL(points = list(free = free, fixed = fixed), candi = candi, schedule = schedule, plotit = TRUE)
objPPL(res)
countPPL(res)
plot(res)
