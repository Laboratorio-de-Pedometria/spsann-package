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
set.seed(2001)
res <- optimMKV(points = 100, candi = candi, covars = covars, maxdist = 500,
                eqn = z ~ dist, vgm = vgm)
objSPSANN(res) # 11.9878
objMKV(points = res, candi = candi, covars = covars, eqn = z ~ dist, vgm = vgm, 
       maxdist = 500)

# 1) GREEDY ALGORITHM ##########################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
vgm <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
set.seed(2001)
res <- optimMKV(points = 100, candi = candi, covars = covars, vgm = vgm, 
                eqn = z ~ dist, greedy = TRUE)
objSPSANN(res) # 11.65344
objMKV(points = res, candi = candi, covars = covars, eqn = z ~ dist, 
       vgm = vgm)

# 2) MANY COVARIATES ###########################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
vgm <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
set.seed(2001)
res <- optimMKV(points = 100, candi = candi, covars = covars, vgm = vgm,
                eqn = z ~ dist + soil + ffreq + x + y)
objSPSANN(res) # 12.03658
objMKV(points = res, candi = candi, covars = covars, 
       eqn = z ~ dist + soil + ffreq + x + y, vgm = vgm)
