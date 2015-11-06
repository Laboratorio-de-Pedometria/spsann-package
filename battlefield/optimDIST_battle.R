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
schedule <- list(initial.acceptance = 0.90, 
                 initial.temperature = 0.5,
                 temperature.decrease = 0.95, chains = 500, 
                 chain.length = 1, stopping = 10)
set.seed(2001)
# \dontrun{
# This example takes more than 5 seconds to run!
res <- optimDIST(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, schedule = schedule, plotit = TRUE)
plot(attr(res, "energy.state"), type = "l")

objSPSANN(res) # 1.6505
objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)
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
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, greedy = TRUE,
                 use.coords = TRUE)
objSPSANN(res) # 1.491634
objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 2) FACTOR COVARIATES WITH THE COORDINATES ####################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 6:7]
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE)
objSPSANN(res) # 1.17612
objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 3) FACTOR AND NUMERIC COVARIATES WITH THE COORDINATES ########################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, c(1, 2, 5:7)]

# Plotting
set.seed(2001)
resA <- optimDIST(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, plotit = TRUE, track = TRUE)
objSPSANN(resA) # 2.924718
objDIST(points = resA, candi = candi, covars = covars, use.coords = TRUE)

# No plotting, but traking
set.seed(2001)
resB <- optimDIST(points = 100, candi = candi, covars = covars, 
                  use.coords = TRUE, track = TRUE)
objSPSANN(resB) # 2.924718
objDIST(points = resB, candi = candi, covars = covars, use.coords = TRUE)

# No plotting, no traking
set.seed(2001)
resC <- optimDIST(points = 100, candi = candi, covars = covars, 
                  use.coords = TRUE)
objSPSANN(resC) # 2.924718
objDIST(points = resC, candi = candi, covars = covars, use.coords = TRUE)
sapply(list(resA, resB, resC), attr, "run")
rm(resA, resB, resC)
gc()
