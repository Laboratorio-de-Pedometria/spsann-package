# Initial settings
rm(list = ls())
gc()
require(pedometrics)
require(gstat)
require(sp)
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
# 0) DEFAULT EXAMPLE ###########################################################
require(sp)
require(gstat)
data(meuse.grid)
candi <- meuse.grid[, 1:2]
covars <- as.data.frame(meuse.grid)
model <- vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
set.seed(2001)
res <- optimMKV(points = 100, candi = candi, covars = covars, maxdist = 500,
                equation = z ~ dist, model = model, iterations = 100,
                plotit = FALSE, track = FALSE, verbose = FALSE)
tail(attr(res, "energy"), 1) # 11.9878
objMKV(points = res, candi = candi, covars = covars, equation = z ~ dist, 
       model = model, maxdist = 500)

# 1) GREEDY ALGORITHM ##########################################################
rm(list = ls())
gc()
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
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
sapply(list.files("src", full.names = TRUE, pattern = ".cpp$"), Rcpp::sourceCpp)
sapply(list.files("R", full.names = TRUE, pattern = ".R$"), source)
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
