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
res <- optimCLHS(points = 100, candi = candi, covars = covars,
                 use.coords = TRUE, iter = 100)
objSPSANN(res) # 
objCLHS(points = res, candi = candi, covars = covars, use.coords = TRUE)
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

# 1) FACTOR COVARIATES USING THE COORDINATES ###################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 6:7]
set.seed(2001)
res <- optimCLHS(points = 100, candi = candi, covars = covars, use.coords = T)
objSPSANN(res) # 
objCLHS(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 2) FACTOR COVARIATES USING THE COORDINATES WITH A FEW POINTS #################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 6:7]
set.seed(2001)
res <-  optimCLHS(points = 10, candi = candi, covars = covars, 
                  use.coords = TRUE, iterations = 200)
objSPSANN(res) # 
objCLHS(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 3) CATEGORICAL COVARIATES WITH MANY COVARIATES AND MANY POINTS ###############
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, rep(c(6, 7), 10)]
set.seed(2001)
res <- optimCLHS(points = 500, candi = candi, covars = covars, use.coords = T)
objSPSANN(res) # 
objCLHS(points = res, candi = candi, covars = covars, use.coords = TRUE)
