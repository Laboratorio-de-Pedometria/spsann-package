# Prepare output
# 
# SUMMARY
# 1. Close progress bar;
# 2. Get last system configuration and its attributes;
# 3. Prepare data on energy states;
# 4. Prepare other attributes (iterations, running time);
# 5. Prepare MOOP attributes (weights, nadir, and utopia);
# 6. Print optimization info and return output.
#
# NOTES
# 1. 
# 2. 
# 
# Progress bar
if (progress) close(pb)

# Last system configuration
res <- new_conf
a <- attributes(res)

# Energy states
if (!track) energies <- new_energy
if (MOOP) {
  criterion <- rbind(energy0, energies)
} else {
  criterion <- c(energy0, energies)
}
a$energy.state <- criterion

# Other attributes
a$iterations <- k
a$running.time <- round(c((proc.time() - time0) / 60)[3], 2)

# MOOP
if (MOOP) {
  a$weights <- weights
  a$nadir <- nadir
  a$utopia <- utopia
}

# Add attributes
attributes(res) <- a

# Print the number of iterations and running time
cat("iterations = ", a$iterations, "\n", sep = "")
cat("running time = ", a$running.time, " minutes", sep = "")

# Output
return (res)
#
# COMMAND
# # Prepare output ###########################################################
# eval(.prepare_output())
# ############################################################################
