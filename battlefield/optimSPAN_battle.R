# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(SpatialTools)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
# PREPARE DATA #################################################################
require(sp)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
nadir <- list(sim = 10, seeds = 1:10)
utopia <- list(user = list(DIST = 0, CORR = 0, PPL = 0, MSSD = 0))
covars <- meuse.grid[, 5]
set.seed(2001)
res <- optimSPAN(points = 100, candi = candi, covars = covars, nadir = nadir,
                 use.coords = TRUE, iterations = 1000, utopia = utopia, 
                 verbose = FALSE, plotit = FALSE, track = FALSE)
tail(attr(res, "energy")$obj, 1) # 
str(res)

# 1) FIRST TEST ################################################################
x11()
set.seed(2001)
tmp <- optimPAN(points = points, candidates = candidates, x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, iterations = iterations, acceptance = acceptance, stopping = stopping, PPL = PPL, ACDC = ACDC, PAN = PAN, plotit = plotit, boundary = boundary, progress = progress, verbose = verbose)
