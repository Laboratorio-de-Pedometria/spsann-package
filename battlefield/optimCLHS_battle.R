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
weights <- list(O1 = 0.5, O3 = 0.5)
schedule <- scheduleSPSANN(chains = 1, initial.temperature = 10)
set.seed(2001)
res <- optimCLHS(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, weights = weights, schedule = schedule)
objSPSANN(res) -
  objCLHS(points = res, candi = candi, covars = covars, use.coords = TRUE, 
          weights = weights)

# 1) FACTOR COVARIATES USING THE COORDINATES ###################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 6:7]
schedule <- scheduleSPSANN(chains = 1, initial.temperature = 10)
set.seed(2001)
res <- optimCLHS(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, schedule = schedule)
objSPSANN(res) -
  objCLHS(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 2) FACTOR COVARIATES USING THE COORDINATES WITH A FEW POINTS #################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 6:7]
schedule <- scheduleSPSANN(chains = 1, initial.temperature = 10)
set.seed(2001)
res <-  optimCLHS(points = 10, candi = candi, covars = covars, 
                  use.coords = TRUE, schedule = schedule, plotit = TRUE)
objSPSANN(res) -
  objCLHS(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 3) CATEGORICAL COVARIATES WITH MANY COVARIATES AND MANY POINTS ###############
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, rep(c(6, 7), 10)]
schedule <- scheduleSPSANN(chains = 1, initial.temperature = 10)
set.seed(2001)
res <- optimCLHS(points = 500, candi = candi, covars = covars, use.coords = T,
                 schedule = schedule, plotit = TRUE)
objSPSANN(res) -
  objCLHS(points = res, candi = candi, covars = covars, use.coords = TRUE)
