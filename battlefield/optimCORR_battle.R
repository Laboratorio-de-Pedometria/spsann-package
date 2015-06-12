# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
# 0) DEFAULT EXAMPLE ###########################################################
require(sp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 5]
set.seed(2001)
res <- optimCORR(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, iterations = 100, plotit = FALSE, 
                 track = FALSE, verbose = FALSE)
tail(attr(res, "energy"), 1) # 0.06386069
objCORR(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 1) FACTOR COVARIATES WITH THE COORDINATES ####################################
rm(list = ls())
gc()
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 6:7]
set.seed(2001)
res <- optimCORR(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, iterations = 100)
tail(attr(res, "energy"), 1) # 2.919182
objCORR(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 2) FACTOR COVARIATES #########################################################
# Tue 9 Jun: objCORR() does not return the same criterion value if 
#            'iterations = 100'
rm(list = ls())
gc()
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, rep(6:7, 3)]
set.seed(2001)
res <- optimCORR(points = 100, candi = candi, covars = covars, iterations = 200)
tail(attr(res, "energy"), 1) # 0.001887334
objCORR(points = res, candi = candi, covars = covars)
