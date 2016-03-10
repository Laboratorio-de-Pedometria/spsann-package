# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(gstat)
require(sp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)

# 0) DEFAULT EXAMPLE ###########################################################
require(sp)
require(gstat)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
vgm <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
schedule <- scheduleSPSANN(initial.temperature = 10, chains = 1,
                           x.max = 1540, y.max = 2060, x.min = 0, 
                           y.min = 0, cellsize = 40)
set.seed(2001)
res <- optimMKV(points = 10, candi = candi, covars = covars, 
                eqn = z ~ dist, vgm = vgm, schedule = schedule)
objSPSANN(res) -
  objMKV(points = res, candi = candi, covars = covars, 
         eqn = z ~ dist, vgm = vgm)

# 1) GREEDY ALGORITHM #########################################################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
vgm <- gstat::vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
schedule <- scheduleSPSANN(initial.temperature = 0.001, initial.acceptance = 0, chains = 1)
set.seed(2001)
res <- optimMKV(
  points = 100, candi = candi, covars = covars, vgm = vgm, eqn = z ~ dist, schedule = schedule, plotit = TRUE, 
  maxdist = 500)
objSPSANN(res) - 
  objMKV(points = res, candi = candi, covars = covars, eqn = z ~ dist, vgm = vgm, maxdist = 500)

# 2) MANY COVARIATES ##########################################################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
vgm <- gstat::vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
schedule <- scheduleSPSANN(initial.temperature = 10, chains = 10)
set.seed(2001)
res <- optimMKV(
  points = 100, candi = candi, covars = covars, vgm = vgm, eqn = z ~ dist + soil + ffreq, plotit = TRUE, 
  schedule = schedule, nmax = 50)
objSPSANN(res) - 
  objMKV(points = res, candi = candi, covars = covars, eqn = z ~ dist + soil + ffreq, vgm = vgm, nmax = 50)

# 3) ORDINARY KRIGING #########################################################################################
# Close to the end of the optimization, the algorithm restarts with the previously best configuration.
# Perhaps this is the reason why the objective value returned by objSPSANN() is not equal to that computed
# with objMKV().
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid, package = "sp")
candi <- meuse.grid[, 1:2]
vgm <- gstat::vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
schedule <- scheduleSPSANN(chains = 500, initial.temperature = 5)
set.seed(2001)
res <- optimMKV(points = 10, candi = candi, vgm = vgm, nmax = 50, plotit = TRUE, schedule = schedule)
objSPSANN(res) - objMKV(points = res, candi = candi, vgm = vgm, nmax = 50)
