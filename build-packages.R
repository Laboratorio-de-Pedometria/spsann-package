################################################################################
# Build packages
################################################################################

# Load packages
library(devtools)

# Build Rcpp files
library(Rcpp)
setwd("~/PROJECTS/pedometrics/pedometrics/pkg/pedometrics") ## isric
setwd("~/PROJECTS/r-packages/pedometrics/pkg/pedometrics") ## laptop
Rcpp::compileAttributes()

# Generate documentation using Rd2oxygen2 ######################################
require(Rd2roxygen)
setwd("~/PROJECTS/pedometrics/pedometrics/pkg/pedometrics") ## isric
setwd("~/PROJECTS/r-packages/pedometrics/pkg/pedometrics") ## laptop
roxygen2::roxygenise()


# Built package
setwd("~/PROJECTS/pedometrics/pedometrics/pkg") ## isric
setwd("~/PROJECTS/r-packages/pedometrics/pkg") ## laptop
system("R CMD build pedometrics")

# Check package
system("R CMD check --help")
system("R CMD check --as-cran pedometrics_0.4-2.tar.gz")



