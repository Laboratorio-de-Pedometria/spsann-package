# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(SpatialTools)
require(magrittr)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)

# 0) DEFAULT EXAMPLE ##########################################################################################
data(meuse.grid, package = 'sp')
candi <- meuse.grid[, 1:2]
schedule <- scheduleSPSANN(
  chains = 1, initial.temperature = 5000000,
  x.max = 1540, y.max = 2060, x.min = 0, y.min = 0, cellsize = 40)
set.seed(2001)
res <- optimMSSD(points = 10, candi = candi, schedule = schedule)
data.frame(
  expected = 247204.8,
  objSPSANN = objSPSANN(res),
  objMSSD = objMSSD(candi = candi, points = res)
)

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
schedule <- scheduleSPSANN(initial.acceptance = 0, initial.temperature = 0.01)
set.seed(2001)
res <- optimMSSD(points = 30, candi = candi, schedule = schedule, plotit = TRUE, boundary = boundary)
data.frame(
  expected = 27982.05,
  objSPSANN = objSPSANN(res),
  objMSSD = objMSSD(candi = candi, points = res)
)

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
objSPSANN(res) - objMSSD(candi = candi, points = res)
plot(res, boundary = boundary)

# 3) USE A COARSER GRID (COMPARED TO CANDI) TO COMPUTE THE OBJECTIVE FUNCITON VALUE ###########################
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
  chains = 100, initial.temperature = 500000, x.max = 1540, y.max = 2060, x.min = 0, y.min = 0, cellsize = 40)
set.seed(2001)
res_candi <- optimMSSD(
  points = 30, candi = candi, schedule = schedule, plotit = TRUE, boundary = boundary)
eval.grid <- sp::spsample(x = boundary, n = 1000, type = 'regular') %>% sp::coordinates()
set.seed(2001)
res_eval_grid <- optimMSSD(
  points = 30, candi = candi, eval.grid = eval.grid, schedule = schedule, plotit = TRUE, boundary = boundary)
objSPSANN(res_candi) - objMSSD(candi = candi, points = res_candi)
objSPSANN(res_eval_grid) - objMSSD(eval.grid = eval.grid, points = res_eval_grid)
objSPSANN(res_candi)
objSPSANN(res_eval_grid)

# 4) THIN AN EXISTING SAMPLE CONFIGURATION ####################################################################
rm(list = ls())
gc()
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid, package = "sp")
boundary <- meuse.grid
sp::coordinates(boundary) <- c("x", "y")
sp::gridded(boundary) <- TRUE
boundary <- rgeos::gUnaryUnion(as(boundary, "SpatialPolygons"))
schedule <- scheduleSPSANN(
  chains = 500, initial.temperature = 6000, x.max = 1540, y.max = 2060, x.min = 0, y.min = 0, cellsize = 0)
eval.grid <- meuse.grid[, 1:2] %>% as.matrix()
candi <- sp::spsample(x = boundary, n = 100, type = 'random') %>% sp::coordinates()
set.seed(2001)
res <- optimMSSD(
  points = 50, candi = candi, eval.grid = eval.grid, schedule = schedule, plotit = TRUE, boundary = boundary)
objSPSANN(res) - objMSSD(points = res, eval.grid = eval.grid)
plot(boundary)
points(res$points[, c("x", "y")], pch = 20)
points(candi, col = 'red')
