# Check suggests
# 
# COMMAND
# eval(.check_suggests())
# 
# SUMMARY
# 1. Check if the suggested packages are installed;
# 2. Compose a message with the suggested packages that are not installed;
# 3. Stop the function call if any suggested package is not installed.
#
# NOTES
# 1. The function requires an object called 'pkg', a vector with the suggested
#    packages.
#    
# Check if suggested packages are installed
id <- !sapply(pkg, requireNamespace, quietly = TRUE)
res <- NULL
if (any(id)) {
  pkg <- paste(pkg[which(id)], collapse = " ")
  res <- paste("Package(s) needed for this function to work but not", "installed: ", pkg, sep = "")
}
if (!is.null(res)) stop (res, call. = FALSE)
