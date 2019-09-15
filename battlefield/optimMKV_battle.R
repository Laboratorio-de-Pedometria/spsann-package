# Initial settings
library(magrittr)
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)

# 0) DEFAULT EXAMPLE ###########################################################
data(meuse.grid, package = "sp")
candi <- meuse.grid[1:1000, 1:2]
covars <- as.data.frame(meuse.grid)[1:1000, ]
vgm <- gstat::vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
schedule <- scheduleSPSANN(
  initial.temperature = 10, chains = 1, x.max = 1540, y.max = 2060, 
  x.min = 0,  y.min = 0, cellsize = 40)
set.seed(2001)
res <- optimMKV(
  points = 10, candi = candi, covars = covars, eqn = z ~ dist, 
  vgm = vgm, schedule = schedule)
objSPSANN(res) - 
  objMKV(points = res, candi = candi, covars = covars,  eqn = z ~ dist, vgm = vgm)

# 1) GREEDY ALGORITHM WITH TOO SMALL NEIGHBOURHOOD SIZE (500 M) ###############################################
# skipped 'singular matrix' error in 'krige'-function
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid, package = 'sp')
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
vgm <- gstat::vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
schedule <- scheduleSPSANN(initial.temperature = 0.001, initial.acceptance = 0, chains = 1, cellsize = 40)
set.seed(2001)
res <- optimMKV(
  points = 100, candi = candi, covars = covars, vgm = vgm, eqn = z ~ dist, schedule = schedule, plotit = TRUE, 
  maxdist = 500)
objSPSANN(res) - 
  objMKV(points = res, candi = candi, covars = covars, eqn = z ~ dist, vgm = vgm, maxdist = 500)

# 2) GREEDY ALGORITHM WITH NEIGHBOURHOOD SET USING NMAX #######################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid, package = 'sp')
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
vgm <- gstat::vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
schedule <- scheduleSPSANN(initial.temperature = 0.001, initial.acceptance = 0, chains = 1, cellsize = 40)
set.seed(2001)
res <- optimMKV(
  points = 100, candi = candi, covars = covars, vgm = vgm, eqn = z ~ dist, schedule = schedule, plotit = TRUE, 
  nmax = 5)
data.frame(
  expected = 13.06498,
  objSPSANN = objSPSANN(res),
  objMKV = objMKV(points = res, candi = candi, covars = covars, eqn = z ~ dist, vgm = vgm, nmax = 5)
)

# 3) MANY COVARIATES ##########################################################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid, package = "sp")
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
vgm <- gstat::vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
schedule <- scheduleSPSANN(initial.temperature = 10, chains = 1, cellsize = 40)
set.seed(2001)
res <- optimMKV(
  points = 100, candi = candi, covars = covars, vgm = vgm, eqn = z ~ dist + soil + ffreq, plotit = TRUE, 
  schedule = schedule, nmax = 50)
data.frame(
  expected = -3.142771e+12,
  objSPSANN = objSPSANN(res),
  objMKV = objMKV(
    points = res, candi = candi, covars = covars, eqn = z ~ dist + soil + ffreq, vgm = vgm, nmax = 50) 
)

# 4) ORDINARY KRIGING #########################################################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid, package = "sp")
candi <- meuse.grid[, 1:2]
vgm <- gstat::vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
schedule <- scheduleSPSANN(chains = 500, initial.temperature = 5, cellsize = 40)
set.seed(2001)
res <- optimMKV(points = 10, candi = candi, eqn = z ~ 1, vgm = vgm, nmax = 50, plotit = TRUE, 
                schedule = schedule)
data.frame(
  expected = 15.9901,
  objSPSANN = objSPSANN(res),
  objMKV = objMKV(points = res, eqn = z ~ 1, candi = candi, vgm = vgm, nmax = 50)
)

# 5) ADD TEN POINTS TO AN EXISTING SAMPLE CONFIGURATION #######################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid, package = "sp")
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid[, c("soil", "dist")])
vgm <- gstat::vgm(psill = 10, model = "Exp", range = 100, nugget = 5)
schedule <- scheduleSPSANN(initial.temperature = 10^200, chains = 10, cellsize = 40)
free <- 10
set.seed(1984)
id <- sample(1:nrow(candi), 90)
fixed <- cbind(id, candi[id, ])
set.seed(2001)
res <- optimMKV(
  points = list(free = free, fixed = fixed), candi = candi, covars = covars, vgm = vgm, 
  eqn = z ~ soil, plotit = TRUE, schedule = schedule)
data.frame(
  fixed = objMKV(points = fixed, candi = candi, covars = covars, eqn = z ~ soil, vgm = vgm),
  expected = 13.62302,
  objSPSANN = objSPSANN(res),
  objMKV = objMKV(points = res, candi = candi, covars = covars, eqn = z ~ soil, vgm = vgm)
)
