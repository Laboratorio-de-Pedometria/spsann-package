# Set progress
#
# COMMAND
# eval(.set_progress())
#
# SUMMARY
# 1. Check if the progress should be monitored.
# 2. Compute the length of the progress bar based on the number and length of
#    the chains, and on the number of points.
# 3. Check if the progress bar should be of type text or tk.
# 4. If a tk progress bar is required, check if the tcltk-package is installed.
# 5. Start counting the processing time.
#
# NOTES
# 1. a tk progress bar is useful when running the algorithm in parallel
#    processors
#
if (!is.null(progress)) {
  max <- n_pts * schedule$chains * schedule$chain.length
  if (progress == "txt") {
    pb <- utils::txtProgressBar(min = 1, max = max, style = 3)
  } else {
    if (progress == "tk") {

      # Check suggests
      pkg <- c("tcltk")
      eval(.check_suggests())

      label <- paste(objective, " with ", n_pts, " points", sep = "")
      pb <- tcltk::tkProgressBar(label = label, min = 1, max = max)
    }
  }
}
# if (progress) {
# pb <- utils::txtProgressBar(min = 1, max = max, style = 3)
# }
time0 <- proc.time()
