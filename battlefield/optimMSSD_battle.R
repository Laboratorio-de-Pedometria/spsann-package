# Initial settings
require(pedometrics)
require(SpatialTools)
require(sp)
source('R/spSANNtools.R')
source('R/spJitter.R')
source('R/optimMSSD.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
Rcpp::sourceCpp('src/calcMSSDCpp.cpp')
Rcpp::sourceCpp('src/updateMSSDCpp.cpp')
# 0) DEFAULT EXAMPLE ###########################################################
require(pedometrics)
require(sp)
require(SpatialTools)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
set.seed(2001)
res <- optimMSSD(points = 100, candi = candi, iterations = 100)
tail(attr(res, "energy.state"), 1) # 11531.03
objMSSD(candi = candi, points = res)

# 1) GREEDY ALGORITHM ##########################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
set.seed(2001)
res <- optimMSSD(points = 100, candi = candi, iterations = 100, greedy = TRUE)
tail(attr(res, "energy.state"), 1) # 13566.74
objMSSD(candi = candi, points = res)
