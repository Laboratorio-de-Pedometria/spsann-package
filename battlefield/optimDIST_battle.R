# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(sp)
source('R/optimDIST.R')
source('R/optimACDC.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
# 0) DEFAULT EXAMPLE ###########################################################
require(pedometrics)
require(sp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 5]
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, iterations = 100)
tail(attr(res, "energy"), 1) # 1.6505
objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 1) GREEDY ALGORITHM ##########################################################
rm(list = ls())
gc()
source('R/optimDIST.R')
source('R/optimACDC.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 5]
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, greedy = TRUE,
                 use.coords = TRUE, iterations = 100)
tail(attr(res, "energy"), 1) # 1.491634
objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 2) FACTOR COVARIATES WITH THE COORDINATES ####################################
rm(list = ls())
gc()
source('R/optimDIST.R')
source('R/optimACDC.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 6:7]
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, iterations = 100)
tail(attr(res, "energy"), 1) # 1.17612
objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 3) FACTOR AND NUMERIC COVARIATES WITH THE COORDINATES ########################
rm(list = ls())
gc()
source('R/optimDIST.R')
source('R/optimACDC.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, c(1, 2, 5:7)]
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, iterations = 100)
tail(attr(res, "energy"), 1) # 2.924718
objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)
