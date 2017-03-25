# Prepare points and candi
# 
# COMMAND
# eval(.prepare_points())
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

# (Fixed) Points ####
if (is(points, "list") && length(points) == 2) {
  points <- points$free
  fixed_pts <- points$fixed
  if (is(fixed_pts, "OptimizedSampleConfiguration")) { # Optimized sample comfiguration
    fixed_pts <- fixed_pts$points
  } 
  if (is.matrix(fixed_pts) || is.data.frame(fixed_pts)) { # Data frame or matrix
    fixed_pts <- as.matrix(fixed_pts)
  } else {
    if (is.integer(fixed_pts) || pedometrics::isNumint(fixed_pts)) {
      if (length(fixed_pts) > 1) { # Integer vector
        fixed_pts <- candi[fixed_pts, ]
      }
      if (length(fixed_pts) == 1) { # Integer value
        stop ("invalid value passed to argument 'points$fixed': integer value")
      }
    }
  }
  n_fixed_pts <- nrow(fixed_pts)
  # Check if 'fixed_pts' has a colunm "id" with the row indexes of 'candi' that correspond to each point
  if (ncol(fixed_pts) != 3 || colnames(fixed_pts)[1] != "id") {
    stop ("missing 'id' column in object 'poinst$fixed'")
  }
}

# (Free) Points ####
if (is(points, "OptimizedSampleConfiguration")) points <- points$points
if (is.matrix(points) || is.data.frame(points)) { # Data frame or matrix
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
# (Fixed) Points ####
if (exists("fixed_pts")) {
points <- rbind(points, fixed_pts)  
}
old_conf <- points
# old_conf <- conf0
#if (COST) cm <- cost[points[, 1]]
