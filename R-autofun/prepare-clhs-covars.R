# Prepare CLHS covariates
#
# COMMAND
# # Prepare 'covars' and create the starting sample matrix 'sm'
# eval(.prepare_clhs_covars())
# 
# SUMMARY
# 1. Check if the coordinates should be used;
# 2. Count the number of covariates;
# 3. Identify the type of covariates (factor, numeric, or both);
# 4. Prepare the starting sample matrix;
# 5. For numeric covariates, compute the break points using continuous sample
#    quantiles, and the population correlation matrix;
# 6. For factor covariates, compute the proportion of population points per
#    marginal factor level.
#
# NOTES
# 1. 

# Use coordinates?
if (use.coords) { covars <- data.frame(covars, candi[, 2:3]) }
n_cov <- ncol(covars)

# Factor and/or numeric?
if (pedometrics::anyFactor(covars)) {
  if (pedometrics::allFactor(covars)) {
    id_fac <- 1:n_cov
    covars_type <- "factor"
    id_num <- NA
  } else {
    id_fac <- which(sapply(covars, is.factor))
    id_num <- which(sapply(covars, is.numeric))
    covars_type <- "both"
  }
} else {
  id_num <- 1:n_cov
  covars_type <- "numeric"
  id_fac <- NA
}

# Sample matrix
sm <- covars[points[, 1], ]

# Break points and population correlation matrix
if (any(covars_type == c("numeric", "both"))) {
  probs <- seq(0, 1, length.out = n_pts + 1)
  breaks <- lapply(covars[, id_num], stats::quantile, probs, na.rm = TRUE)
  pcm <- stats::cor(x = covars[, id_num], use = "complete.obs")
}

# Proportion of population points per marginal factor level
if (any(covars_type == c("factor", "both"))) {
  pop_prop <- lapply(covars[, id_fac], function(x) table(x) / n_candi)
}
