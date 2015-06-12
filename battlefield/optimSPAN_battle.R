# Initial settings
rm(list = ls())
gc()
require(ASRtools)
require(pedometrics)
require(sp)
require(rgeos)
require(Hmisc)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
# PREPARE DATA #################################################################
require(sp)
data(meuse.grid)
candidates <- meuse.grid[, 1:2]
coordinates(candidates) <- ~ x + y
gridded(candidates) <- TRUE
boundary <- as(candidates, "SpatialPolygons")
boundary <- gUnionCascaded(boundary)
candidates <- coordinates(candidates)
candidates <- matrix(cbind(c(1:dim(candidates)[1]), candidates), ncol = 3)
points <- 100
x.max <- diff(bbox(boundary)[1, ])
y.min <- x.min <- 40
y.max <- diff(bbox(boundary)[2, ])
iterations <- 10000
acceptance <- list(initial = 0.99, cooling = iterations / 10)
weights <- list(strata = 0.5, correl = 0.5)
stopping <- list(max.count = iterations / 10)
plotit <- TRUE
progress <- TRUE
verbose <- TRUE
PPL <- list(lags.type  = "exponential",
            criterion  = "distribution",
            lags       = 7,
            lags.base  = 2,
            pre.distri = NULL,
            cutoff     = sqrt((x.max * x.max) + (y.max * y.max)))
ACDC <- list(covars      = meuse.grid[, 1],
             covars.type = "numeric",
             use.coords  = TRUE,
             strata.type = "equal.area",
             weights     = list(strata = 0.5, correl = 0.5))
PAN <- list(weights = list(PPL = 1/3, ACDC = 1/3, MSSD = 1/3),
            nadir   = list(sim = 10, save.sim = FALSE, user = NULL, 
                           abs = NULL))
# 1) FIRST TEST ################################################################
x11()
set.seed(2001)
tmp <- optimPAN(points = points, candidates = candidates, x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, iterations = iterations, acceptance = acceptance, stopping = stopping, PPL = PPL, ACDC = ACDC, PAN = PAN, plotit = plotit, boundary = boundary, progress = progress, verbose = verbose)
