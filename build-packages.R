################################################################################
# Build packages
################################################################################

# Load packages
library(devtools)

# Build Rcpp files
library(Rcpp)
Rcpp::compileAttributes()

# Generate documentation using Rd2oxygen2 ######################################
require(Rd2roxygen)
roxygen2::roxygenise()


# Built and check package
setwd("~/PROJECTS/r-packages")
system("R CMD build spsann")
system("R CMD check --as-cran spsann_0.0.0.9000.tar.gz")
