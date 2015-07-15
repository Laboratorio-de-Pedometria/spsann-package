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
covars <- meuse.grid[, 5]
set.seed(2001)
# \dontrun{
# This example takes more than 5 seconds to run!
res <- optimCORR(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE)
objSPSANN(res) # 0.06386069
objCORR(points = res, candi = candi, covars = covars, use.coords = TRUE)
# }
# Random sample
pts <- sample(1:nrow(candi), 5)
pts <- cbind(pts, candi[pts, ])
objCORR(points = pts, candi = candi, covars = covars, use.coords = TRUE)
# 1) FACTOR COVARIATES WITH THE COORDINATES ####################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 6:7]
set.seed(2001)
res <- optimCORR(points = 100, candi = candi, covars = covars, use.coords = TRUE)
objSPSANN(res) # 2.919182
objCORR(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 2) FACTOR COVARIATES #########################################################
# Tue 9 Jun: objCORR() does not return the same criterion value if 
#            'iterations = 100'
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, rep(6:7, 3)]
set.seed(2001)
res <- optimCORR(points = 100, candi = candi, covars = covars, iterations = 200)
objSPSANN(res) # 0.001887334
objCORR(points = res, candi = candi, covars = covars)
