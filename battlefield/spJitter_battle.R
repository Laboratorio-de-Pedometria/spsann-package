rm(list = ls())
gc()
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
# Default example ##############################################################
require(sp)
data(meuse.grid)
meuse.grid <- as.matrix(meuse.grid[, 1:2])
meuse.grid <- matrix(cbind(1:dim(meuse.grid)[1], meuse.grid), ncol = 3)
pts1 <- sample(c(1:dim(meuse.grid)[1]), 155)
pts2 <- meuse.grid[pts1, ]
pts3 <- spJitterFinite(points = pts2, candi = meuse.grid, x.min = 40,
                       x.max = 100, y.min = 40, y.max = 100, which.point = 10)
plot(meuse.grid[, 2:3], asp = 1, pch = 15, col = "gray")
points(pts2[, 2:3], col = "red", cex = 0.5)
points(pts3[, 2:3], pch = 19, col = "blue", cex = 0.5)
