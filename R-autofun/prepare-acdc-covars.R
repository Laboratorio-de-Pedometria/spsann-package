# Prepare ACDC covariates and base data
#
# SUMMARY
# 1. Check if the coordinates should be used;
# 2. Identify the type of covariates (factor or numeric);
# 3. If needed, convert numeric covariates to factor covariates;
# 4. Count the number of covariates;
# 5. Prepare the starting sample matrix.
#
# NOTES
# 1. 

# Use coordinates?
if (use.coords) {
  covars <- data.frame(covars, candi[, 2:3])
}

# Factor or numeric?
covars.type <- ifelse(pedometrics::anyFactor(covars), "factor", "numeric")

# Convert numeric covariates to factor covariates
if (covars.type == "factor") {
  if (!pedometrics::allFactor(covars)) {
    i <- which(sapply(covars, is.factor) == FALSE)
    mes <- paste("converting ", length(i), 
                 " numeric covariates to factor covariates", sep = "")
    message(mes)
    num_covars <- data.frame(covars[, i])
    breaks <- .strataACDC(n.pts = n_pts, strata.type = strata.type, 
                          covars = num_covars, covars.type = "numeric")[[1]]
    num_covars <- pedometrics::cont2cat(x = num_covars, breaks = breaks)
    covars[, i] <- num_covars
  }
}
n_cov <- ncol(covars)
sm <- covars[points[, 1], ]

# COMMAND
# # Prepare 'covars' and create the starting sample matrix 'sm'
# eval(.prepare_acdc_covars())
