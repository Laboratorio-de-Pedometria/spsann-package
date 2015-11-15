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
criterion <- rbind(energy0, energies)
# if (MOOP) {
#   criterion <- rbind(energy0, energies)
# } else {
#   criterion <- c(energy0, energies)
# }
a$energy.state <- criterion

# Other attributes
a$iterations <- k

# Running time
rt <- as.numeric(c(proc.time() - time0)[3])
if (rt > 3600) {
  rt <- list(time = round(rt / 3600, 2), unit = "hours")
} else {
  if (rt > 60) {
    rt <- list(time = round(rt / 60, 2), unit = "minutes")
  } else {
    rt <- list(time = round(rt, 2), unit = "seconds")
  }
}
a$running.time <- rt

# MOOP
if (MOOP) {
  a$weights <- weights
  if (objective != "CLHS") {
    a$nadir <- nadir
    a$utopia <- utopia 
  }
}

# Add attributes
attributes(res) <- a

# Print the number of iterations and running time
cat("iterations = ", a$iterations, "\n", sep = "")
cat("running time = ", rt$time, " ", rt$unit, sep = "")

# Output
return (res)
#
# COMMAND
# # Prepare output
# eval(.prepare_output())
