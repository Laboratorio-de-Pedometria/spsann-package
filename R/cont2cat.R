#' Continuous to categorical
#' 
#' Convert a continuous variable into a categorical one.
#' 
#' @param x Data frame or matrix with the continuous variables to be converted
#' into categorical variables (factors).
#' 
#' @param breaks List with the lower and upper limits to be used to break the
#' continuous variable. Using a list allows breaking the variables into a
#' different number of classes.
#' 
#' @return
#' data.frame
#' 
# MAIN FUNCTION ################################################################
cont2cat <-
  function (x, breaks) {
    for (i in 1:ncol(x)) {
      x[, i] <- cut2(x[, i], breaks[[i]])
    }
    return (x)
  }
