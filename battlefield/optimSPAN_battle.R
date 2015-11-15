# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(SpatialTools)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)

# 0) DEFAULT EXAMPLE ###########################################################
require(sp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
nadir <- list(sim = 10, seeds = 1:10)
utopia <- list(user = list(DIST = 0, CORR = 0, PPL = 0, MSSD = 0))
covars <- meuse.grid[, 5]
schedule <- scheduleSPSANN(chains = 1, initial.temperature = 1)
set.seed(2001)
res <- optimSPAN(points = 100, candi = candi, covars = covars, nadir = nadir,
                 use.coords = TRUE, utopia = utopia, schedule = schedule)
objSPSANN(res) -
  objSPAN(points = res@points, candi = candi, covars = covars, nadir = nadir,
          use.coords = TRUE, utopia = utopia)
