# Generated by autofun (0.0.0.9000): do not edit by hand!!!
# Please edit source code in spsann-package/R-autofun/prepare-output.R
.prepare_output<-function(...){
expression(if (!is.null(progress)) close(pb), rt <- as.numeric(c(proc.time() - time0)[3]), 
    if (rt > 3600) {
      rt <- list(time = round(rt / 3600, 2), unit = "hours")
    } else {
      if (rt > 60) {
        rt <- list(time = round(rt / 60, 2), unit = "minutes")
      } else {
        rt <- list(time = round(rt, 2), unit = "seconds")
      }
    }, if (!track) energies <- new_energy, energies <- rbind(energy0, energies), 
    res <- list(
      
      # The optimized sample configuration.
      # The content of this data frame MUST be numeric.
      points = data.frame(new_conf, free = c(rep(1, n_pts), rep(0, n_fixed_pts))),
      
      # Information about the spatial simulated annealing
      spsann = list(
        acceptance = data.frame(initial = NA_real_),
        cellsize = data.frame(x = NA_real_, y = NA_real_),
        chains = data.frame(total = NA_integer_, used = NA_integer_, length = NA_integer_),
        jitter = data.frame(x = rep(NA_real_, 2), y = rep(NA_real_, 2), row.names = c("min", "max")),
        running = data.frame(time = NA_real_, units = NA_character_), 
        stopping = NA_integer_,
        temperature = data.frame(initial = NA_real_, final = NA_real_)),
      
      # Information about the objective function
      objective = list(
        name = NA_character_,
        energy = data.frame(NA_real_),
        nadir = data.frame(NA_real_),
        utopia = data.frame(NA_real_),
        weights = data.frame(NA_real_))
    ), class(res) <- "OptimizedSampleConfiguration", res$spsann <- list(
      acceptance = data.frame(initial = schedule$initial.acceptance),
      cellsize = data.frame(x = cellsize[1], y = cellsize[2]),
      chains = data.frame(total = schedule$chains, used = i, length = schedule$chain.length),
      jitter = data.frame(x = c(x.min, x_max0), y = c(y.min, y_max0), row.names = c("min", "max")),
      running = data.frame(time = rt[[1]], units = rt[[2]]),
      stopping = schedule$stopping,
      temperature = data.frame(initial = schedule$initial.temperature, final = actual_temp)
    ), res$objective <- list(
      name = objective,
      energy = energies,
      nadir = if (MOOP && objective != "CLHS") data.frame(nadir),
      utopia = if (MOOP && objective != "CLHS") data.frame(utopia),
      weights = if (MOOP) data.frame(weights)
    ), if (objective %in% c("ACDC", "CLHS", "CORR", "DIST", "SPAN")) {
      if (objective != "CLHS") res$objective$strata.type <- strata.type
      res$objective$use.coords <- use.coords
    }, if (objective %in% c("PPL", "SPAN")) {
      res$objective$lags <- lags
      res$objective$lags.type <- lags.type
      res$objective$lags.base <- lags.base
      # res$objective$cutoff <- cutoff
      res$objective$cutoff <- tail(lags, 1)
      res$objective$criterion <- criterion
      res$objective$distri <- distri # delete
      res$objective$target.distribution <- distri
      res$objective$final.distribution <- ppl
      res$objective$pairs <- pairs
    }, if (objective == "MKV") {
      res$objective$eqn <- eqn
      res$objective$vgm <- vgm
      res$objective$krige.stat <- krige.stat
    }, cat("running time = ", rt$time, " ", rt$unit, sep = ""), 
    return (res))
}

