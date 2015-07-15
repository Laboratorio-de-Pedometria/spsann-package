# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(SpatialTools)
require(sp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
# 0) DEFAULT EXAMPLE ###########################################################
require(sp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
set.seed(2001)
# \dontrun{
# This example takes more than 5 seconds to run!
res <- optimMSSD(points = 100, candi = candi)
objSPSANN(res) # 11531.03
objMSSD(candi = candi, points = res)
# }
# Random sample
pts <- sample(1:nrow(candi), 5)
pts <- cbind(pts, candi[pts, ])
objMSSD(candi = candi, points = pts)
# 1) GREEDY ALGORITHM ##########################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
set.seed(2001)
res <- optimMSSD(points = 100, candi = candi, iterations = 100, greedy = TRUE)
objSPSANN(res) # 13566.74
objMSSD(candi = candi, points = res)
