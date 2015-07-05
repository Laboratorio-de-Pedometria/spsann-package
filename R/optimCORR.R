#' Optimization of sample configurations for spatial trend identification and
#' estimation
#'
#' Optimize a sample configuration for spatial trend identification and 
#' estimation. A criterion is defined so that the sample reproduces the 
#' bivariate association/correlation between the covariates (\bold{CORR}).
#'
#' @template spJitter_doc
#' @template spSANN_doc
#' @template ACDC_doc
#' @template MOOP_doc
#' 
#' @details
#' The \emph{correlation} between two numeric covariates is measured using the 
#' Pearson's r, a descriptive statistic that ranges from -1 to +1. 
#' This statistic is also known as the linear correlation coefficient.
#' 
#' When the set of covariates includes factor covariates, any numeric covariate 
#' is transformed into a factor covariate. The factor levels are defined 
#' using the marginal sampling strata created using one of the two methods 
#' available (equal-area or equal-range strata) (see more details at
#' \code{optimDist()}).
#' 
#' The \emph{association} between two factor covariates is measured using the 
#' Cramér's v, a descriptive statistic that ranges from 0 to +1. The closer to 
#' +1 the Cramér's v is, the stronger the association between two factor 
#' covariates. The main weakness of using the Cramér's v is that, while the 
#' Pearson's r shows the degree and direction of the association between two 
#' covariates (negative or positive), the Cramér's v only measures the degree 
#' of association (weak or strong).
#'
#' @return
#' \code{optimCORR} returns a matrix: the optimized sample configuration.
#' 
#' \code{objCORR} returns a numeric value: the energy state of the sample
#' configuration - the objective function value.
#'
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[clhs]{clhs}}, \code{\link[pedometrics]{cramer}}
#' @importFrom pedometrics cramer
#' @importFrom pedometrics is.numint
#' @importFrom pedometrics cont2cat
#' @importFrom SpatialTools dist2
#' @aliases optimCORR objCORR
#' @export
#' @examples
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' covars <- meuse.grid[, 5]
#' set.seed(2001)
#' res <- optimCORR(points = 100, candi = candi, covars = covars, 
#'                  use.coords = TRUE, iterations = 100, plotit = FALSE, 
#'                  track = FALSE, verbose = FALSE)
#' tail(attr(res, "energy"), 1) # 0.06386069
#' objCORR(points = res, candi = candi, covars = covars, use.coords = TRUE)
# MAIN FUNCTION ################################################################
optimCORR <-
  function (
    # CORR
    covars, use.coords = FALSE, strata.type = "area", 
    # SPSANN
    points, candi, x.max, x.min, y.max, y.min, iterations,
    acceptance = list(initial = 0.99, cooling = iterations / 10),
    stopping = list(max.count = iterations / 10), plotit = TRUE,
    boundary, progress = TRUE, verbose = TRUE, greedy = FALSE, track = TRUE, 
    # MOOP
    weights = NULL, nadir = NULL, utopia = NULL) {
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Check other arguments
    check <- .optimACDCcheck(candi = candi, covars = covars, 
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Prepare 'covars' and create the starting sample matrix 'sm'
    eval(.prepare_acdc_covars())
    
    # Base data and initial energy state
    pcm <- .corCORR(obj = covars, covars.type = covars.type)
    scm <- .corCORR(obj = sm, covars.type = covars.type)
    energy0 <- .objCORR(scm = scm, pcm = pcm)

    # Other settings for the simulated annealing algorithm
    old_scm <- scm
    best_scm <- scm
    old_sm <- sm
    new_sm <- sm
    best_sm <- sm
    count <- 0
    old_energy <- energy0
    best_energy <- Inf
    if (progress) pb <- txtProgressBar(min = 1, max = iterations, style = 3)
    time0 <- proc.time()

    # Begin the main loop
    for (k in 1:iterations) {
      
      # Plotting and jittering
      eval(.plot_and_jitter())
      
      # Update sample and correlation matrices, and energy state
      new_sm[wp, ] <- covars[new_conf[wp, 1], ]
      new_scm <- .corCORR(obj = new_sm, covars.type = covars.type)
      new_energy <- .objCORR(scm = new_scm, pcm = pcm)
      
      # Evaluate the new system configuration
      if (greedy) {
        random_prob <- 1
      } else {
        random_prob <- runif(1)
      }
      actual_prob <- acceptance[[1]] * exp(-k / acceptance[[2]])
      if (track) accept_probs[k] <- actual_prob
      if (new_energy <= old_energy) {
        old_conf <- new_conf
        old_energy <- new_energy
        count <- 0
        old_sm <- new_sm
        old_scm <- new_scm
      } else {
        if (new_energy > old_energy & random_prob <= actual_prob) {
          old_conf <- new_conf
          old_energy <- new_energy
          count <- count + 1
          old_sm <- new_sm
          old_scm <- new_scm
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... p = ",
                random_prob, "\n")
          }
        } else {
          new_energy <- old_energy
          new_conf <- old_conf
          count <- count + 1
          new_sm <- old_sm
          new_scm <- old_scm
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... stops at",
                stopping[[1]], "\n")
          }
        }
      }
      # Best energy state
      if (track) energies[k] <- new_energy
      if (new_energy < best_energy / 1.0000001) {
        best_k <- k
        best_conf <- new_conf
        best_energy <- new_energy
        best_old_energy <- old_energy
        old_conf <- old_conf
        best_sm <- new_sm
        best_old_sm <- old_sm
        best_scm <- new_scm
        best_old_scm <- old_scm
      }
      
      # Freezing parameters
      if (count == stopping[[1]]) {
        if (new_energy > best_energy * 1.000001) {
          old_conf <- old_conf
          new_conf <- best_conf
          old_energy <- best_old_energy
          new_energy <- best_energy
          count <- 0
          new_sm <- best_sm
          new_scm <- best_scm
          old_sm <- best_old_sm
          old_scm <- best_old_scm
          cat("\n", "reached maximum count with suboptimal configuration\n")
          cat("\n", "restarting with previously best configuration\n")
          cat("\n", count, "iteration(s) with no improvement... stops at",
              stopping[[1]], "\n")
        } else {
          break
        }
      }
      if (progress) setTxtProgressBar(pb, k)
    }
    
    # Prepare output
    eval(.prepare_output())
  }
# FUNCTION - CALCULATE ENERGY STATE ############################################
#' @rdname optimCORR
#' @export
objCORR <-
  function (points, candi, covars, use.coords = FALSE, strata.type = "area") {
    
    # Check other arguments
    check <- .optimACDCcheck(candi = candi, covars = covars, 
                             use.coords = use.coords, strata.type = strata.type)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare 'covars' and create the starting sample matrix 'sm'
    eval(.prepare_acdc_covars())
    
    # Calculate the energy state
    pcm <- .corCORR(obj = covars, covars.type = covars.type)
    scm <- .corCORR(obj = sm, covars.type = covars.type)
    energy <- .objCORR(scm = scm, pcm = pcm)
    
    return (energy)
  }
# INTERNAL FUNCTION - CALCULATE ASSOCIATION/CORRELATION MATRIX ################
# obj: the matrix of covariates, for the population correlation matrix, and the
#      sample matrix, for the sample correlation matrix
# covars.type: the type of covariate
.corCORR <-
  function (obj, covars.type) {
    
    if (covars.type == "numeric") { 
      cor(x = obj, use = "complete.obs")
      
    } else { # Factor covariates
      pedometrics::cramer(x = obj)
    }
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE ############################
# scm: sample correlation matrix
# pcm: population correlation matrix
.objCORR <-
  function (scm, pcm) {
    sum(abs(pcm - scm))
  }
