# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(sp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)

# 0) DEFAULT EXAMPLE ###########################################################
# \dontrun{
# This example takes more than 5 seconds to run!
require(sp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
schedule <- scheduleSPSANN(chains = 1, initial.temperature = 30)
set.seed(2001)
res <- optimPPL(points = 100, candi = candi, schedule = schedule)
objSPSANN(res) - objPPL(points = res, candi = candi)
countPPL(points = res, candi = candi)
# }

# 1) Point pairs with many chains ##############################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
schedule <- scheduleSPSANN(chains = 500, initial.temperature = 500)
set.seed(2001)
res <- optimPPL(points = 100, candi = candi, pairs = TRUE, schedule = schedule,
                plotit = TRUE)
objSPSANN(res) - objPPL(points = res, pairs = TRUE, candi = candi)
countPPL(points = res, candi = candi, pairs = TRUE)

# 2) Points per lag - select sample points from candi ##########################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]

# random selection
points <- 100
set.seed(2001)
res <- countPPL(points = points, candi = candi, cutoff = 1000)
set.seed(2001)
objPPL(points = points, candi = candi, cutoff = 1000)
sum(points - res$ppl) # 346
# vector of indexes
points <- 1:100
res <- countPPL(points = points, candi = candi, cutoff = 1000)
objPPL(points = points, candi = candi, cutoff = 1000)
sum(length(points) - res$ppl) # 266

# 3) Unit test #################################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
set.seed(2001)
countPPL(points = 100, candi = candi, lags = 1, cutoff = Inf)
set.seed(2001)
countPPL(points = 100, candi = candi, lags = 1, pairs = TRUE, cutoff = Inf)
