# 0) DEFAULT EXAMPLE ##########################################################################################
rm(list = ls())
gc()
devtools::load_all()
# Example
data(meuse.grid, package = "sp")
candi <- meuse.grid[, 1:2]
nadir <- list(sim = 10, seeds = 1:10)
utopia <- list(user = list(DIST = 0, CORR = 0, PPL = 0, MSSD = 0))
covars <- meuse.grid[, 5]
schedule <- scheduleSPSANN(
  chains = 1, initial.acceptance = c(0.9, 0.99),
  initial.temperature = 0.85, x.max = 1540, y.max = 2060, 
  x.min = 0, y.min = 0, cellsize = 40)
weights <- list(CORR = 1/6, DIST = 1/6, PPL = 1/3, MSSD = 1/3)
set.seed(2001, sample.kind = "Round")
res <- optimSPAN(
  points = 10, candi = candi, covars = covars, nadir = nadir, weights = weights, 
  use.coords = TRUE, utopia = utopia, schedule = schedule)
objSPSANN(res) - objSPAN(
  points = res, candi = candi, covars = covars, nadir = nadir, weights = weights,
  use.coords = TRUE, utopia = utopia)

# 1) ADD TEN POINTS TO EXISTING SAMPLE CONFIGURATION ##########################################################
rm(list = ls())
gc()
devtools::load_all()
data(meuse.grid, package = "sp")
candi <- meuse.grid[, 1:2]
nadir <- list(sim = 10, seeds = 1:10)
utopia <- list(user = list(DIST = 0, CORR = 0, PPL = 0, MSSD = 0))
covars <- meuse.grid[, 5]
schedule <- scheduleSPSANN(
  chains = 100, initial.acceptance = c(0.9, 0.99),
  initial.temperature = 0.1, x.max = 1540, y.max = 2060, 
  x.min = 0, y.min = 0, cellsize = 40)
free <- 10
set.seed(1984, sample.kind = "Round")
id <- sample(1:nrow(candi), 40)
fixed <- cbind(id, candi[id, ])
objSPAN(
  points = fixed, candi = candi, covars = covars, nadir = nadir, use.coords = TRUE, utopia = utopia,
  weights = list(CORR = 1/6, DIST = 1/6, PPL = 1/3, MSSD = 1/3))
set.seed(2001)
res <- optimSPAN(
  points = list(free = free, fixed = fixed), candi = candi, covars = covars, nadir = nadir,
  use.coords = TRUE, utopia = utopia, schedule = schedule, plotit = TRUE,
  weights = list(CORR = 1/6, DIST = 1/6, PPL = 1/3, MSSD = 1/3))
objSPSANN(res) - 
  objSPAN(
  points = res, candi = candi, covars = covars, nadir = nadir, use.coords = TRUE, utopia = utopia,
  weights = rev(list(CORR = 1/6, DIST = 1/6, PPL = 1/3, MSSD = 1/3)))
