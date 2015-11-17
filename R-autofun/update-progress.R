# Update progress bar
#
# COMMAND
# eval(.update_progress())
# 
# SUMMARY
# 1. Check if there is a progress bar to be updated
# 2. Check if the progress bar is of type txt or tk
#
# NOTES
# 1. tk progress bar should be useful when running the algorithm in parallel
#    processors
# 
if (!is.null(progress)) {
  if (progress == "txt") {
    utils::setTxtProgressBar(pb, k)
  } else {
    if (progress == "tk") {
      tcltk::setTkProgressBar(pb, k)
    }
  }
}
# if (progress) {
# utils::setTxtProgressBar(pb, k)
# }
