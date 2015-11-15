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
schedule <- scheduleSPSANN(initial.temperature = 0.5, chains = 1)
set.seed(2001)
# \dontrun{
# This example takes more than 5 seconds to run!
res <- optimCORR(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, schedule = schedule)
objSPSANN(res) -
  objCORR(points = res@points, candi = candi, covars = covars, 
          use.coords = TRUE)
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
schedule <- scheduleSPSANN(initial.temperature = 0.5, chains = 10)
set.seed(2001)
res <- optimCORR(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, schedule = schedule, plotit = TRUE)
objSPSANN(res) # 3.744234
objCORR(points = res@points, candi = candi, covars = covars, use.coords = TRUE)
