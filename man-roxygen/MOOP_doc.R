#  Template documentation for multi-objective optimization problems (MOOP)
####################################################################################################
#' @param nadir List with named sub-arguments. Three options are available:
#'  * `sim`: the number of simulations that should be used to estimate the nadir point, and
#'    `seeds` vector defining the random seeds for each simulation;
#'  * `user`: a list of user-defined nadir values named after the respective objective
#'    functions to which they apply;
#'  * `abs`: logical for calculating the nadir point internally (experimental).
#'
#' @param weights List with named sub-arguments. The weights assigned to each one of the objective
#' functions that form the multi-objective combinatorial optimization problem. They must be named
#' after the respective objective function to which they apply. The weights must be equal to or
#' larger than 0 and sum to 1.
#'
#' @param utopia List with named sub-arguments. Two options are available:
#'  * `user`: a list of user-defined values named after the respective objective functions to which
#'  they apply;
#'  * `abs`: logical for calculating the utopia point internally (experimental).
#'
#' @details
#' The help page of [spsann::minmaxPareto()] contains details on how __spsann__ solves the
#' multi-objective combinatorial optimization problem of finding a globally optimum sample
#' configuration that meets multiple, possibly conflicting, sampling objectives.
