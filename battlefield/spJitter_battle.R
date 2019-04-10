# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(sp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)

# Default example #############################################################################################
require(sp)
data(meuse.grid)
meuse.grid <- as.matrix(meuse.grid[, 1:2])
meuse.grid <- matrix(cbind(1:dim(meuse.grid)[1], meuse.grid), ncol = 3)
pts1 <- sample(c(1:dim(meuse.grid)[1]), 155)
pts2 <- meuse.grid[pts1, ]
pts3 <- spJitter(points = pts2, candi = meuse.grid, x.min = 40,
                 x.max = 100, y.min = 40, y.max = 100, 
                 which.point = 10, cellsize = 40)
plot(meuse.grid[, 2:3], asp = 1, pch = 15, col = "gray")
points(pts2[, 2:3], col = "red", cex = 0.5)
points(pts3[, 2:3], pch = 19, col = "blue", cex = 0.5)

# Cluster of points
pts1 <- c(1:55)
pts2 <- meuse.grid[pts1, ]
pts3 <- spJitter(points = pts2, candi = meuse.grid, x.min = 40,
                 x.max = 80, y.min = 40, y.max = 80, 
                 which.point = 1, cellsize = 40)
plot(meuse.grid[, 2:3], asp = 1, pch = 15, col = "gray")
points(pts2[, 2:3], col = "red", cex = 0.5)
points(pts3[, 2:3], pch = 19, col = "blue", cex = 0.5)

# TWO POINTS, SHORT DISTANCE ##################################################################################
rm(list = ls())
gc()
require(pedometrics)
require(sp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
meuse.grid <- as.matrix(meuse.grid[1:96, 1:2])
meuse.grid <- matrix(cbind(1:dim(meuse.grid)[1], meuse.grid), ncol = 3)
pts1 <- sample(c(1:dim(meuse.grid)[1]), 2)
pts2 <- meuse.grid[pts1, ]
pts3 <- spJitter(
  points = pts2, candi = meuse.grid, x.min = 0, x.max = 40, y.min = 0, y.max = 40, which.point = 1, 
  cellsize = 40)
plot(meuse.grid[, 2:3], asp = 1, pch = 15, col = "gray")
points(pts2[, 2:3], col = "red")
points(pts3[, 2:3], pch = 19, col = "blue", cex = 0.5)
pts2[, 2:3] - pts3[, 2:3]

# FINITE SET OF CANDIDATE LOCATIONS (cellsize = 0) ############################################################
# Should not have duplicated sample points
rm(list = ls())
gc()
require(pedometrics)
require(sp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
data(meuse.grid)
schedule <- spsann::scheduleSPSANN(
  chains = 5, initial.temperature = 10, temperature.decrease = 0.5, x.max = 1540, y.max = 2060, x.min = 0, 
  y.min = 0, cellsize = 0)
set.seed(2001)
res <- spsann::optimCLHS(
  points = 400, candi = meuse.grid[, c('x', 'y')], covars = meuse.grid[, c(1, 5)], 
  use.coords = F, schedule = schedule, track = F, plotit = T, weights = list(O1 = 1/2, O3 = 1/2))
res.pt <- res$points
sp::coordinates(res.pt) <- ~ x + y
sp::zerodist(res.pt)
