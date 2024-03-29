# Generated by autofun (0.0.0.9000): do not edit by hand!!!
# Please edit source code in spsann-package/R-autofun/prepare-clhs-covars.R
.prepare_clhs_covars<-function(...){
expression(if (use.coords) { covars <- data.frame(covars, candi[, 2:3]) }, 
    n_cov <- ncol(covars), if (pedometrics::anyFactor(covars)) {
      if (pedometrics::allFactor(covars)) {
        id_fac <- 1:n_cov
        covars_type <- "factor"
        id_num <- NA
      } else {
        id_fac <- which(sapply(covars, is.factor))
        id_num <- which(sapply(covars, is.numeric))
        covars_type <- "both"
      }
    } else {
      id_num <- 1:n_cov
      covars_type <- "numeric"
      id_fac <- NA
    }, sm <- covars[points[, 1], ], if (any(covars_type == c("numeric", "both"))) {
      probs <- seq(0, 1, length.out = n_pts + n_fixed_pts + 1)
      breaks <- lapply(covars[, id_num], stats::quantile, probs, na.rm = TRUE)
      pcm <- stats::cor(x = covars[, id_num], use = "complete.obs")
    }, if (any(covars_type == c("factor", "both"))) {
      # pop_prop <- lapply(covars[, id_fac], function(x) table(x) / n_candi)
      pop_count <- lapply(covars[, id_fac], function(x) table(x))
    })
}

