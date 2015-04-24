################################################################################
# Build packages
################################################################################

# Load packages
library(devtools)
library(Rcpp)
require(Rd2roxygen)

# Build Rcpp files
Rcpp::compileAttributes()

# Generate documentation using Rd2oxygen2 ######################################
roxygen2::roxygenise()

# Built and check package ######################################################
setwd("~/PROJECTS/r-packages") # Laptop
setwd("~/alessandro") # ISRIC desktop

system("R CMD build spsann")
system("R CMD check --as-cran spsann_0.0.0.9004.tar.gz")

setwd("~/alessandro/spsann") # ISRIC desktop
