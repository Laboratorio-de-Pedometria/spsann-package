# Prepare points and candi
#
# SUMMARY
# 1. Compute 'n_candi', the number of finite candidate locations;
# 2. Add the 'id' column to 'candi, and force 'candi' to be of class matrix;
# 3. Prepare the 'points' for the spatial simulated annealing routine;
# 4. Compute 'n_pts', the number of points;
# 5. Define 'conf0' and 'old_conf', objects that are used later in the 
#    optimization;
# 6. FUTURE FEATURE: Create 'cm', the starting cost matrix (if needed).
#
# NOTES
# 1. This code chunk is used in both 'optim...()' and 'obj...()' functions.
#    The family of 'obj...()' functions may not use the argument 'candi'. This
#    is the reason for evaluating if 'candi' is missing.
if (!missing(candi)) {
  n_candi <- nrow(candi)
  candi <- as.matrix(cbind(id = 1:n_candi, candi))
}

# Points
if (is(points, "OptimizedSampleConfiguration")) points <- points@points
if (is.matrix(points) || is.data.frame(points)) { # Data frame of matrix
  points <- as.matrix(points)
} else {
  if (is.integer(points) || pedometrics::isNumint(points)) {
    if (length(points) > 1) { # Integer vector
      points <- candi[points, ]
    }
    if (length(points) == 1) { # Integer value
      points <- sample(1:n_candi, points)
      points <- candi[points, ] 
    }
  }
}
n_pts <- nrow(points)
conf0 <- points
old_conf <- conf0
#if (COST) cm <- cost[points[, 1]]
#
# COMMAND
# # Prepare points and candi
# eval(.prepare_points())
