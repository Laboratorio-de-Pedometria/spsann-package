# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(sp)
source('R/optimMKV.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
# 0) DEFAULT EXAMPLE ###########################################################
require(pedometrics)
require(sp)
require(gstat)
require(plyr)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
model <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
set.seed(2001)
res <- optimMKV(points = 100, candi = candi, covars = covars, 
                equation = z ~ dist, model = model, iterations = 100)
tail(attr(res, "energy"), 1) # 11.61896
objMKV(points = res, candi = candi, covars = covars, equation = z ~ dist, 
       model = model)

# 1) GREEDY ALGORITHM ##########################################################
rm(list = ls())
gc()
source('R/optimMKV.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
model <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
set.seed(2001)
res <- optimMKV(points = 100, candi = candi, covars = covars, model = model, 
                equation = z ~ dist, iterations = 100, greedy = TRUE)
tail(attr(res, "energy"), 1) # 11.65344
objMKV(points = res, candi = candi, covars = covars, equation = z ~ dist, 
       model = model)

# 2) MANY COVARIATES ###########################################################
rm(list = ls())
gc()
source('R/optimMKV.R')
source('R/spSANNtools.R')
source('R/spJitter.R')
Rcpp::sourceCpp('src/spJitterCpp.cpp')
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
model <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
set.seed(2001)
res <- optimMKV(points = 100, candi = candi, covars = covars, model = model,
                equation = z ~ dist + soil + ffreq + x + y, plotit = FALSE,
                iterations = 100)
tail(attr(res, "energy"), 1) # 12.03658
objMKV(points = res, candi = candi, covars = covars, 
       equation = z ~ dist + soil + ffreq + x + y, model = model)
