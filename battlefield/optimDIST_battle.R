# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(sp)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
# 0) DEFAULT EXAMPLE ###########################################################
require(sp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, 5]
set.seed(2001)
res <- optimDIST(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, iterations = 100, plotit = FALSE, 
                 track = FALSE, verbose = FALSE)
tail(attr(res, "energy"), 1) # 1.6505
objDIST(points = res, candi = candi, covars = covars, use.coords = TRUE)

# 1) GREEDY ALGORITHM ##########################################################
rm(list = ls())
gc()
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
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
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
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
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- meuse.grid[, c(1, 2, 5:7)]

# Plotting
set.seed(2001)
resA <- optimDIST(points = 100, candi = candi, covars = covars, 
                 use.coords = TRUE, iterations = 100, verbose = FALSE)
tail(attr(resA, "energy"), 1) # 2.924718
objDIST(points = resA, candi = candi, covars = covars, use.coords = TRUE)

# No plotting, but traking
set.seed(2001)
resB <- optimDIST(points = 100, candi = candi, covars = covars, 
                  use.coords = TRUE, iterations = 100, verbose = FALSE,
                  plotit = FALSE)
tail(attr(resB, "energy"), 1) # 2.924718
objDIST(points = resB, candi = candi, covars = covars, use.coords = TRUE)

# No plotting, no traking
set.seed(2001)
resC <- optimDIST(points = 100, candi = candi, covars = covars, 
                  use.coords = TRUE, iterations = 100, verbose = FALSE,
                  plotit = FALSE, track = FALSE)
tail(attr(resC, "energy"), 1) # 2.924718
objDIST(points = resC, candi = candi, covars = covars, use.coords = TRUE)

sapply(list(resA, resB, resC), attr, "run")
rm(resA, resB, resC)
gc()
