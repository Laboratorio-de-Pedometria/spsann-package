# Prepare ACDC covariates
#
# SUMMARY
# 1. Identify the type of covariates (factor or numeric);
# 2. Prepare the object 'covars';
# 2.1 Check if the coordinates should be used;
# 2.2 Convert numeric covariates to factor covariates if needed;
# 3. Count the number of covariates;
# 4. Prepare the starting sample matrix.
#
# NOTES
# 1. 
covars.type <- ifelse(pedometrics::is.any.factor(covars), "factor", "numeric")
# covars <- .covarsACDC(covars = covars, covars.type = covars.type, 
#                       candi = candi, use.coords = use.coords, n.pts = n_pts, 
#                       strata.type = strata.type)
# Use coordinates?
if (use.coords) {
  covars <- data.frame(covars, candi[, 2:3])
}
if (covars.type == "factor") {
  
  # Convert numeric covariates to factor covariates
  if (!pedometrics::is.all.factor(covars)) {
    i <- which(sapply(covars, is.factor) == FALSE)
    mes <- paste("converting ", length(i), 
                 " numeric covariates to factor covariates", sep = "")
    message(mes)
    num_covars <- data.frame(covars[, i])
    breaks <- .numStrata(n.pts = n_pts, covars = num_covars, 
                         strata.type = strata.type)[[1]]
    num_covars <- pedometrics::cont2cat(x = num_covars, breaks = breaks)
    covars[, i] <- num_covars
  }
}
n_cov <- ncol(covars)
sm <- covars[points[, 1], ]

# COMMAND
# # Prepare 'covars' and create the starting sample matrix 'sm' ##############
# eval(.prepare_acdc_covars())
# ############################################################################
