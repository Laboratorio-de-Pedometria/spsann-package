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
schedule <- scheduleSPSANN(chains = 1, initial.temperature = 5000)
set.seed(2001)
# \dontrun{
# This example takes more than 5 seconds to run!
res <- optimMSSD(points = 100, candi = candi, schedule = schedule)
objSPSANN(res) - objMSSD(candi = candi, points = res)
# }
# Random sample
pts <- sample(1:nrow(candi), 5)
pts <- cbind(pts, candi[pts, ])
objMSSD(candi = candi, points = pts)
# 1) GREEDY ALGORITHM WITH MANY CHAINS #########################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
schedule <- scheduleSPSANN(chains = 500, initial.acceptance = 0, 
                           initial.temperature = 0.01)
set.seed(2001)
res <- optimMSSD(points = 100, candi = candi, schedule = schedule, 
                 plotit = TRUE)
objSPSANN(res) - objMSSD(candi = candi, points = res)
