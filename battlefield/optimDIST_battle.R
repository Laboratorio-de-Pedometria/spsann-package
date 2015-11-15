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
res <- optimDIST(points = 100, candi = candi, covars = covars,
                 use.coords = TRUE, schedule = schedule)
objSPSANN(res) -
  objDIST(points = res@points, candi = candi, covars = covars, 
          use.coords = TRUE)
# }
# Random sample
pts <- sample(1:nrow(candi), 5)
pts <- cbind(pts, candi[pts, ])
objDIST(points = pts, candi = candi, covars = covars, use.coords = TRUE)

# 1) GREEDY ALGORITHM ##########################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 5]
schedule <- scheduleSPSANN(initial.temperature = 1e-10)
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars,
                 use.coords = TRUE, schedule = schedule, plotit = TRUE)
objSPSANN(res) # 1.492452
objDIST(points = res@points, candi = candi, covars = covars, use.coords = TRUE)

# 2) FACTOR COVARIATES WITH THE COORDINATES ####################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 6:7]
schedule <- scheduleSPSANN(initial.temperature = 0.5, chains = 1)
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, schedule = schedule, plotit = TRUE)
objSPSANN(res) # 1.570828
objDIST(points = res@points, candi = candi, covars = covars, use.coords = TRUE)

# 3) FACTOR AND NUMERIC COVARIATES WITH THE COORDINATES ########################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, c(1, 2, 5:7)]
schedule <- scheduleSPSANN(initial.temperature = 0.5, chains = 1)
set.seed(2001)
resA <- optimDIST(points = 100, candi = candi, covars = covars, 
                  use.coords = TRUE, plotit = TRUE, schedule = schedule)
objSPSANN(resA) # 3.313735
objDIST(points = resA@points, candi = candi, covars = covars, use.coords = TRUE)
