# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(sp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)

# Default example ##############################################################
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

# TWO POINTS, SHORT DISTANCE ###################################################
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
pts3 <- spJitter(points = pts2, candi = meuse.grid, x.min = 0,
                 x.max = 40, y.min = 0, y.max = 40, 
                 which.point = 1, cellsize = 40)
plot(meuse.grid[, 2:3], asp = 1, pch = 15, col = "gray")
points(pts2[, 2:3], col = "red")
points(pts3[, 2:3], pch = 19, col = "blue", cex = 0.5)
pts2[, 2:3] - pts3[, 2:3]
