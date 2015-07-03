# Supress R CMD check note 'no visible binding for global variable ...'
# Source: http://stackoverflow.com/a/12429344/3365410
if (getRversion() >= "2.15.1") {
  utils::globalVariables(names = c("covars.type", "sm", "n_pts", "n_cov", 
                                   "n_candi", "wp", "conf0", "pre.distri",
                                   "pop.prop"))
}
