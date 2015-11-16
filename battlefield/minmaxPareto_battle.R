# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(sp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)

# Default example ##############################################################
require(sp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, c(1, 2)]

# CORR
schedule <- scheduleSPSANN(initial.acceptance = 0.1, chains = 1, 
                           x.max = 1540, y.max = 2060, x.min = 0, 
                           y.min = 0, cellsize = 40)
set.seed(2001)
osc_corr <- optimCORR(points = 10, candi = candi, covars = covars, 
                      schedule = schedule)

# DIST
set.seed(2001)
osc_dist <- optimDIST(points = 10, candi = candi, covars = covars,
                      schedule = schedule)

# PPL
set.seed(2001)
osc_ppl <- optimPPL(points = 10, candi = candi, schedule = schedule)

# MSSD
set.seed(2001)
osc_mssd <- optimMSSD(points = 10, candi = candi, schedule = schedule)

# Pareto
pareto <- minmaxPareto(osc = list(DIST = osc_dist, CORR = osc_corr,
                                  PPL = osc_ppl, MSSD = osc_mssd),
                       candi = candi, covars = covars)
pareto
