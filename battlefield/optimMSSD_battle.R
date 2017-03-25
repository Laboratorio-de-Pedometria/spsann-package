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
schedule <- scheduleSPSANN(chains = 1, initial.temperature = 5000000,
                           x.max = 1540, y.max = 2060, x.min = 0, 
                           y.min = 0, cellsize = 40)
set.seed(2001)
res <- optimMSSD(points = 10, candi = candi, schedule = schedule)
objSPSANN(res) - objMSSD(candi = candi, points = res)

# 1) GREEDY ALGORITHM WITH MANY CHAINS ########################################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid, package = "sp")
boundary <- meuse.grid
sp::coordinates(boundary) <- c("x", "y")
sp::gridded(boundary) <- TRUE
boundary <- rgeos::gUnaryUnion(as(boundary, "SpatialPolygons"))
candi <- meuse.grid[, 1:2]
schedule <- scheduleSPSANN(
  initial.acceptance = 0, initial.temperature = 0.01)
set.seed(2001)
res <- optimMSSD(
  points = 30, candi = candi, schedule = schedule, plotit = TRUE,
  boundary = boundary)
objSPSANN(res)
objMSSD(candi = candi, points = res)
plot(res, boundary = boundary)

# 2) ADD THREE POINTS TO AN ALREADY EXISTING SPATIAL SAMPLE CONFIGURATION #####################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid, package = "sp")
boundary <- meuse.grid
sp::coordinates(boundary) <- c("x", "y")
sp::gridded(boundary) <- TRUE
boundary <- rgeos::gUnaryUnion(as(boundary, "SpatialPolygons"))
candi <- meuse.grid[, 1:2]
schedule <- scheduleSPSANN(
  chains = 50000, initial.temperature = 500000, x.max = 1540, y.max = 2060, x.min = 0, y.min = 0, 
  cellsize = 40)
set.seed(2001)
fixed <- candi[sample(1:nrow(candi), 10), ]
free <- 3
set.seed(1984)
res <- optimMSSD(
  points = list(fixed = fixed, free = free), candi = candi, schedule = schedule, plotit = TRUE, 
  boundary = boundary)
res$points[-(1:free), ]
fixed
objSPSANN(res)
objMSSD(candi = candi, points = res)
