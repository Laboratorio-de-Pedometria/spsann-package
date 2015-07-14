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
nadir <- list(sim = 10, seeds = 1:10)
utopia <- list(user = list(DIST = 0, CORR = 0))
covars <- meuse.grid[, 5]
set.seed(2001)
res <- optimACDC(points = 100, candi = candi, covars = covars, nadir = nadir,
                 use.coords = TRUE, utopia = utopia)
objSPSANN(res) # 0.5272031
objACDC(points = res, candi = candi, covars = covars, use.coords = TRUE, 
        nadir = nadir, utopia = utopia)
# MARGINAL DISTRIBUTION
par(mfrow = c(3, 3))
# Covariates
i <- sample(1:nrow(candi), 100)
hist(candi[, 1], breaks = 10)
hist(candi[, 2], breaks = 10)
hist(covars, breaks = 10)
# Optimized sample
hist(candi[res[, 1], 1], breaks = 10)
hist(candi[res[, 1], 2], breaks = 10)
hist(covars[res[, 1]], breaks = 10)
# Random sample
hist(candi[i, 1], breaks = 10)
hist(candi[i, 2], breaks = 10)
hist(covars[i], breaks = 10)

# LINEAR CORRELATION
# Covariates
cor(cbind(candi[, 1], candi[, 2], covars))
# Optimized sample
cor(cbind(candi[res[, 1], 1], candi[res[, 1], 2], covars[res[, 1]]))
# Random sample
cor(cbind(candi[i, 1], candi[i, 2], covars[i]))

# 1) FACTOR COVARIATES USING THE COORDINATES, WITH USER-DEFINED NADIR ##########
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
nadir <- list(user = list(DIST = 10, CORR = 1))
utopia <- list(user = list(DIST = 0, CORR = 0))
covars <- meuse.grid[, 6:7]
set.seed(2001)
res <- optimACDC(points = 100, candi = candi, covars = covars, nadir = nadir, 
                 use.coords = TRUE, utopia = utopia)
objSPSANN(res) # 1.552125
objACDC(points = res, candi = candi, covars = covars, use.coords = TRUE, 
        nadir = nadir, utopia = utopia)

# 2) FACTOR COVARIATES USING THE COORDINATES WITH A FEW POINTS #################
# Tue 9 Jun: objACDC() does not return the same criterion value if 
#            'iterations = 100'
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 6:7]
nadir <- list(sim = 10, seeds = 1:10)
utopia <- list(user = list(CORR = 0, DIST = 0))
set.seed(2001)
res <- optimACDC(points = 10, candi = candi, covars = covars, nadir = nadir,
                 use.coords = TRUE, iterations = 200, utopia = utopia)
objSPSANN(res) # 0.7908377
objACDC(points = res, candi = candi, covars = covars, use.coords = TRUE,
        nadir = nadir, utopia = utopia)

# 3) CATEGORICAL COVARIATES WITH MANY COVARIATES AND MANY POINTS ###############
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, rep(c(6, 7), 10)]
nadir <- list(sim = 10, seeds = 1:10)
utopia <- list(user = list(CORR = 0, DIST = 0))
set.seed(2001)
res <- optimACDC(points = 500, candi = candi, covars = covars, nadir = nadir, 
                 use.coords = TRUE, utopia = utopia)
objSPSANN(res) # 0.620825
objACDC(points = res, candi = candi, covars = covars, use.coords = TRUE, 
        nadir = nadir, utopia = utopia)

# 4) NUMERIC COVARIATES USING THE COORDINATES, WITH USER-DEFINED NADIR #########
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
nadir <- list(user = list(DIST = 10, CORR = 1))
utopia <- list(user = list(DIST = 0, CORR = 0))
covars = meuse.grid[, 5]
set.seed(2001)
res <- optimACDC(points = 100, candi = candi, covars = covars, nadir = nadir, 
                 use.coords = TRUE, utopia = utopia)
objSPSANN(res) # 0.1340229
objACDC(points = res, candi = candi, covars = covars, nadir = nadir,
        use.coords = TRUE, utopia = utopia)
