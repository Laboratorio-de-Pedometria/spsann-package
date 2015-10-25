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
set.seed(2001)
res <- optimPPL(points = 100, candi = candi)
objSPSANN(res) # 160
objPPL(points = res, candi = candi)
countPPL(points = res, candi = candi)

# 1) Point pairs ###############################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
set.seed(2001)
res <- optimPPL(points = 100, candi = candi, pairs = TRUE)
objSPSANN(res) # 7059.143
objPPL(points = res, pairs = TRUE, candi = candi)
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
res <- countPPL(points = points, candi = candi)
set.seed(2001)
objPPL(points = points, candi = candi)
sum(points - res$ppl) # 216
# vector of indexes
points <- 1:100
res <- countPPL(points = points, candi = candi)
objPPL(points = points, candi = candi)
sum(length(points) - res$ppl) # 200

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

# 4) GREEDY ALGORITH ###########################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
set.seed(2001)
res <- optimPPL(points = 100, candi = candi, greedy = TRUE)
objSPSANN(res) # 165
objPPL(points = res, candi = candi)
countPPL(points = res, candi = candi)

# 5) OLD INFINITE LOOP #########################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
nx <- ny <- 100
x <- 1:nx
x <- x - 0.5
y <- x
candi <- expand.grid(x = x, y = y)
coordinates(candi) <- ~ x + y
gridded(candi) <- TRUE
boundary <- as(candi, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candi <- coordinates(candi)
x.max <- diff(bbox(boundary)[1, ])
y.max <- diff(bbox(boundary)[2, ])
cutoff <- sqrt((x.max * x.max) + (y.max * y.max))# / 2
points <- 49
set.seed(2001)
sample_b <- optimPPL(points = points, candi = candi, cutoff = cutoff, 
                     x.max = x.max, x.min = 2, y.max = y.max, 
                     y.min = 2, boundary = boundary, iterations = 7000, 
                     pairs = TRUE, greedy = TRUE)
objSPSANN(sample_b) # 1176
objPPL(sample_b, cutoff = cutoff, pairs = TRUE)
countPPL(sample_b, cutoff = cutoff, pairs = TRUE)
