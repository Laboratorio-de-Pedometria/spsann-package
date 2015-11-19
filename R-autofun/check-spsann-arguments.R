# Check spsann arguments
# 
# SUMMARY
# 1. Check the mandatory OPTIM arguments 'points' and 'candi';
# 2. If all is fine, check the mandatory OPTIM argument 'iterations';
# 3. If all is fine, check the mandatory OPTIM argument 'acceptance';
# 4. If all is fine, check the mandatory OPTIM argument 'stopping';
# 5. Check if the optional OPTIM arguments 'weights', 'nadir' and 'utopia' are 
#    missing;
# 6. If the optional OPTIM arguments are not missing, check the mandatory MOOP
#    arguments 'weights', 'nadir' and 'utopia';
# 7. Check if acceptance probabilities should be tracked.
#
# NOTES
# 1. This code chunk is used only with the family of 'optim...()' functions.
# 2. A good reading on argument checking is the article at
#    http://www.r-bloggers.com/a-warning-about-warning/
# 
res <- NULL

# Check missing arguments
aa <- c("points", "candi")
bb <- c(missing(points), missing(candi))
if (any(bb)) {
  i <- which(bb == TRUE)
  res <- c("missing argument: '", aa[i], "'\n", sep = "")
} else {
  
  # Argument 'candi'
  aa <- any(apply(candi, 2, is.numeric) == FALSE)
  bb <- ncol(candi) != 2
  cc <- any(c(c("x", "y") != colnames(candi)) == TRUE)
  if (aa || bb || cc) {
    res <- c("'candi' must have two named numeric columns: 'x' and 'y'")
  } else {
    
    # Argument 'schedule'
    if (!is.list(schedule) || length(schedule) != 11) {
      res <- c("'schedule' must be a list with 11 sub-arguments")
    }
  }
}

# Multi-objective optimization
aa <- objective %in% c("ACDC", "CLHS", "SPAN")
# aa <- all(c(!is.null(weights), !is.null(utopia), !is.null(nadir)))
MOOP <- ifelse(aa, TRUE, FALSE)
# if (aa) {
#   MOOP <- length(which(weights != 0))
#   MOOP <- ifelse(MOOP > 1, TRUE, FALSE)
#   #COST <- ifelse(weights$COST == 0, FALSE, TRUE)
# }

# Argument 'weights'
if (MOOP) {
  aa <- !is.list(weights)
  bb <- is.null(names(weights))
  cc <- length(weights) < 2
  if (aa || bb || cc) {
    res <- c("'weights' must be a list with named components")
  } else {
    aa <- sum(unlist(weights)) != 1
    if (aa) {
      res <- c("'weights' must sum to 1")
    }
  }
}

# Arguments 'nadir' and 'utopia'
if (MOOP && objective != "CLHS") {
  aa <- !is.list(utopia)
  bb <- !length(utopia) == 1
  cc <- is.null(names(utopia))
  if (aa || bb || cc) {
    res <- c("'utopia' must be a list with one named component")
  } else {
    
    aa <- !is.list(nadir)
    if (aa) {
      res <- c("'nadir' must be a list with named components")
    }
    aa <- names(nadir)
    if (length(aa) >= 2) {
      if (length(aa) > 2) {
        res <- c("you must choose a single 'nadir' option")
      }
    } else {
      if (aa == "sim" || aa == "seeds") { 
        res <- c("the number of simulations and their seeds must be set")
      }
      if (aa == "abs") {
        res <- c("sorry but 'nadir' cannot be calculated")
      }
    }  
  }
}

# Argument 'track'
if (track || plotit) {
  if (plotit) { track <- TRUE }
  if (MOOP) {
    energies <- as.data.frame(
      matrix(NA_real_, nrow = 1, ncol = length(weights) + 1))
    colnames(energies) <- c("obj", names(weights))
  } else {
    energies <- data.frame(obj = NA_real_)
    # energies <- vector()
  }
}

# Output and cleanup
if (!is.null(res)) stop (res, call. = FALSE)
# rm(aa, bb, cc, dd, res)
rm(aa, bb, cc, res)
#
# COMMAND
# # Check spsann arguments ###################################################
# eval(.check_spsann_arguments())
# ############################################################################
