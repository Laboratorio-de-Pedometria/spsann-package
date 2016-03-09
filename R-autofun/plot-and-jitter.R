# Plotting and jittering
# 
# COMMAND
# eval(.plot_and_jitter())
# 
# SUMMARY
# 1. Select one point to be jittered;
# 2. Plot the evlolution of the energy state, and the old and new system
#    configurations possibly indicating the selected point;
# 3. Jitter the selected point;
# 4. Update x.max and y.max
# 5. ...
#
# NOTES
# 1. 
# 2. 
# 
# if (acceptance$by == "iterations") {
#   wp <- sample(1:n_pts, 1)
# } else {
#   wp <- wp + 1
# }
if (plotit && pedometrics::isNumint(k / 10)) {
  .spSANNplot(
    energy0 = energy0, energies = energies, k = k, boundary = boundary, new_conf = new_conf[, 2:3], 
    # acceptance = acceptance, accept_probs = accept_probs,
    # wp = wp,
    conf0 = conf0[, 2:3], y_max0 = y_max0, y.max = y.max, x_max0 = x_max0, x.max = x.max, 
    best.energy = best_energy, best.k = best_k, MOOP = MOOP)
}
new_conf <- spJitter(
  points = old_conf, candi = candi, x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min,
  #finite = finite,
  cellsize = cellsize, which.point = wp)
# if (acceptance$by == "iterations") {
#   x.max <- x_max0 - (k / iterations) * (x_max0 - x.min)
#   y.max <- y_max0 - (k / iterations) * (y_max0 - y.min)
# } else {
#   if (wp == n_pts) {
#     chain <- chain + 1
#     x.max <- x_max0 - (chain / n_chains) * (x_max0 - x.min)
#     y.max <- y_max0 - (chain / n_chains) * (y_max0 - y.min)
#   }
# }
# if (COST) {
#   new_row <- cost[new_conf[wp, 1]]
#   new_cm[wp] <- new_row
# }
#
