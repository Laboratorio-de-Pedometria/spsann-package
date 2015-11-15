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

# Energy states
if (!track) energies <- new_energy
criterion <- rbind(energy0, energies)

res <- new("OptimizedSampleConfiguration", points = data.frame(new_conf))
slot(res, "spsann") <- 
  list(acceptance = data.frame(initial = schedule$initial.acceptance),
       cellsize = data.frame(x = cellsize[1], y = cellsize[2]),
       chains = data.frame(total = schedule$chains, used = i, 
                           length = schedule$chain.length),
       jitter = data.frame(x = c(x.min, x_max0), y = c(y.min, y_max0), 
                           row.names = c("min", "max")),
       running = data.frame(time = rt[[1]], units = rt[[2]]),
       stopping = schedule$stopping,
       temperature = data.frame(initial = schedule$initial.temperature,
                                final = actual_temp))
slot(res, "objective") <-
  list(name = objective,
       energy = criterion,
       nadir = if (MOOP && objective != "CLHS") data.frame(nadir),
       utopia = if (MOOP && objective != "CLHS") data.frame(utopia),
       weights = if (MOOP) data.frame(weights))

if (objective %in% c("ACDC", "CLHS", "CORR", "DIST", "SPAN")) {
  if (objective != "CLHS") res@objective$strata.type <- strata.type
  res@objective$use.coords <- use.coords
}
if (objective %in% c("PPL", "SPAN")) {
  res@objective$lags <- lags
  res@objective$lags.type <- lags.type
  res@objective$lags.base <- lags.base
  res@objective$cutoff <- cutoff
  res@objective$criterion <- criterion
  res@objective$distri <- distri
  res@objective$pairs <- pairs
}
if (objective == "MKV") {
  res@objective$eqn <- eqn
  res@objective$vgm <- vgm
  res@objective$krige.stat <- krige.stat
}

# Print the running time
cat("running time = ", rt$time, " ", rt$unit, sep = "")

# Output
return (res)
#
# COMMAND
# # Prepare output
# eval(.prepare_output())
