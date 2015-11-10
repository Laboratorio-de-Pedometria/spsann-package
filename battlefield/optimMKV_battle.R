# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(gstat)
require(sp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)

# 0) DEFAULT EXAMPLE ###########################################################
# \dontrun{
# # This example takes more than 5 seconds to run!
require(sp)
require(gstat)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
vgm <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
schedule <- scheduleSPSANN(initial.temperature = 0.5, chains = 1)
set.seed(1984)
res <- optimMKV(points = 100, candi = candi, covars = covars, maxdist = 600,
                eqn = z ~ dist, vgm = vgm, schedule = schedule)
objSPSANN(res) -
  objMKV(points = res, candi = candi, covars = covars, eqn = z ~ dist, 
         vgm = vgm, maxdist = 600)
# }

# 1) GREEDY ALGORITHM ##########################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
vgm <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
schedule <- scheduleSPSANN(initial.temperature = 0.001, initial.acceptance = 0,
                           chains = 1)
set.seed(2001)
res <- optimMKV(points = 100, candi = candi, covars = covars, vgm = vgm, 
                eqn = z ~ dist, schedule = schedule, plotit = TRUE, 
                maxdist = 1000)
objSPSANN(res) # 11.76088
objMKV(points = res, candi = candi, covars = covars, eqn = z ~ dist, 
       vgm = vgm, maxdist = 1000)

# 2) MANY COVARIATES ###########################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
vgm <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
schedule <- scheduleSPSANN(initial.temperature = 0.5, chains = 1)
set.seed(2001)
res <- optimMKV(points = 100, candi = candi, covars = covars, vgm = vgm,
                eqn = z ~ dist + soil + ffreq + x + y, plotit = TRUE)
objSPSANN(res) # 12.05515
objMKV(points = res, candi = candi, covars = covars, 
       eqn = z ~ dist + soil + ffreq + x + y, vgm = vgm)
