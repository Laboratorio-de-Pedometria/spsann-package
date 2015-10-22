#' Optimization of sample configurations for spatial trend identification
#' and estimation (IV)
#'
#' Optimize a sample configuration for spatial trend identification and 
#' estimation using the method proposed by Minasny and McBratney (2006), known 
#' as the conditioned Latin hypercube sampling. An utility function \emph{U} is
#' defined so that the sample reproduces the marginal distribution and
#' correlation matrix of the numeric covariates, and the class proportions of 
#' the factor covariates (\bold{CLHS}). The utility function is obtained
#' aggregating three objective functions: \bold{O1}, \bold{O2}, and \bold{O3}.
#'
#' @template spJitter_doc
#' @template spSANN_doc
#' @inheritParams optimACDC
#' 
#' @section Marginal sampling strata:
#' Reproducing the marginal distribution of the numeric covariates depends upon
#' the definition of marginal sampling strata. \emph{Equal-area} marginal 
#' sampling strata are defined using the sample quantiles estimated with 
#' \code{\link[stats]{quantile}} using a continuous function (\code{type = 7}),
#' that is, a function that interpolates between existing covariate values to 
#' estimate the sample quantiles -- this is the procedure implemented in the 
#' method of Minasny and McBratney (2006), which creates breakpoints that do 
#' not occur in the population of existing covariate values. Depending on the 
#' level of discretization of the covariate values, that is, how many 
#' significant digits they have, this can create repeated breakpoints, 
#' resulting in empty marginal sampling strata. The number of empty marginal 
#' sampling strata will ultimately depend on the frequency distribution of the 
#' covariate, and on the number of sampling points.
#' 
#' @section Correlation between numeric covariates:
#' The \emph{correlation} between two numeric covariates is measured using the 
#' sample Pearson's \emph{r}, a descriptive statistic that ranges from $-1$ to 
#' $+1$. This statistic is also known as the sample linear correlation 
#' coefficient.
#' 
#' @section Multi-objective optimization:
#' A method of solving a multi-objective optimization problem (MOOP) is to 
#' aggregate the objective functions into a single \emph{utility function U}. 
#' In the \pkg{spsann} package, as in the original CLHS, the aggregation is 
#' performed using the \emph{weighted sum method}, which uses weights to 
#' incorporate the preferences of the user about the relative importance of 
#' each objective function. When the user has no preference, the objective 
#' functions receive equal weights.
#' 
#' The weighted sum method is affected by the relative magnitude of the 
#' different objective function values. The objective functions implemented in 
#' \code{optimCLHS} have different units and orders of magnitude. The 
#' consequence is that the objective function with the largest values, 
#' generally \bold{O1}, may have a numerical dominance during the optimization. 
#' In other words, the weights will not express the true preferences of the 
#' user, and the meaning of the utility function becomes unclear -- the 
#' optimization will favour the objective function which is numerically
#' dominant.
#' 
#' An efficient solution to avoid numerical dominance is to transform the 
#' objective functions so that they are constrained to the same approximate 
#' range of values, at least in the end of the optimization. However, as in the 
#' original CLHS, \code{optimCLHS} uses the naive aggregation method, which 
#' ignores that the three objective functions have different units and orders 
#' of magnitude. The same aggregation procedure is implemented in the 
#' \pkg{clhs} package.
#' 
#' @return
#' \code{optimCLHS} returns a matrix: the optimized sample configuration.
#' 
#' \code{objCLHS} returns a numeric value: the energy state of the sample
#' configuration - the objective function value.
#' 
#' @references
#' Minasny, B.; McBratney, A. B. A conditioned Latin hypercube method for
#' sampling in the presence of ancillary information. \emph{Computers &
#' Geosciences}, v. 32, p. 1378-1388, 2006.
#'
#' Minasny, B.; McBratney, A. B. Conditioned Latin Hypercube Sampling for
#' calibrating soil sensor data to soil properties. Chapter 9. Viscarra Rossel,
#' R. A.; McBratney, A. B.; Minasny, B. (Eds.) \emph{Proximal Soil Sensing}.
#' Amsterdam: Springer, p. 111-119, 2010.
#'
#' Roudier, P.; Beaudette, D.; Hewitt, A. A conditioned Latin hypercube sampling
#' algorithm incorporating operational constraints. \emph{5th Global Workshop on
#' Digital Soil Mapping}. Sydney, p. 227-231, 2012.
#'
#' @note
#' The (only) difference of the \code{optimCLHS} function to the original 
#' Fortran implementation of Minasny and McBratney (2006), and to the 
#' \code{clhs} function implemented in the \pkg{\link[clhs]{clhs}} package by
#' Pierre Roudier, is in the annealing schedule.
#'
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[clhs]{clhs}}, \code{\link[spsann]{optimACDC}}
#' @concept spatial trend
#' @examples
#' \dontrun{
#' # This example takes more than 5 seconds to run!
#' require(sp)
#' data(meuse.grid)
#' candi <- meuse.grid[, 1:2]
#' covars <- meuse.grid[, 5]
#' set.seed(2001)
#' res <- optimCLHS(points = 100, candi = candi, covars = covars, 
#'                  use.coords = TRUE, iter = 100)
#' objSPSANN(res) # 
#' objCLHS(points = res, candi = candi, covars = covars, use.coords = TRUE) 
#' 
#' # MARGINAL DISTRIBUTION
#' par(mfrow = c(3, 3))
#' # Covariates
#' i <- sample(1:nrow(candi), 100)
#' hist(candi[, 1], breaks = 10)
#' hist(candi[, 2], breaks = 10)
#' hist(covars, breaks = 10)
#' # Optimized sample
#' hist(candi[res[, 1], 1], breaks = 10)
#' hist(candi[res[, 1], 2], breaks = 10)
#' hist(covars[res[, 1]], breaks = 10)
#' # Random sample
#' hist(candi[i, 1], breaks = 10)
#' hist(candi[i, 2], breaks = 10)
#' hist(covars[i], breaks = 10)
#' 
#' # LINEAR CORRELATION
#' # Covariates
#' cor(cbind(candi[, 1], candi[, 2], covars))
#' # Optimized sample
#' cor(cbind(candi[res[, 1], 1], candi[res[, 1], 2], covars[res[, 1]]))
#' # Random sample
#' cor(cbind(candi[i, 1], candi[i, 2], covars[i]))
#' }
# MAIN FUNCTION ################################################################
optimCLHS <-
  function (points, candi, iterations = 50000, 
    # O1, O2, and O3
    covars, use.coords = FALSE, 
    # SPSANN
    x.max, x.min, y.max, y.min,
    acceptance = list(initial = 0.80, cooling = iterations / 10),
    stopping = list(max.count = iterations / 10), plotit = FALSE, track = FALSE,
    boundary, progress = TRUE, verbose = FALSE, greedy = FALSE,
    # MOOP
    weights = list(O1 = 1/3, O2 = 1/3, O3 = 1/3)) {
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Check other arguments
    check <- 
      .optimCLHScheck(candi = candi, covars = covars, use.coords = use.coords)
    if (!is.null(check)) { stop (check, call. = FALSE) }
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Prepare 'covars' and create the starting sample matrix 'sm'
    # Use coordinates?
    if (use.coords) { covars <- data.frame(covars, candi[, 2:3]) }
    n_cov <- ncol(covars)
    # Factor or numeric?
    if (pedometrics::anyFactor(covars)) {
      if (pedometrics::allFactor(covars)) {
        id_fac <- 1:n_cov
        covars_type <- "factor"
      } else {
        id_fac <- which(sapply(covars, is.factor))
        id_num <- which(sapply(covars, is.numeric))
        covars_type <- "both"
      }
    } else {
      id_num <- 1:n_cov
      covars_type <- "numeric"
    }
    sm <- covars[points[, 1], ]
    
    # Base data and initial energy state
    if (any(covars_type == c("numeric", "both")) {
      
      # O3
      pcm <- stats::cor(x = covars[, id_num], use = "complete.obs")
      scm <- stats::cor(x = sm[, id_num], use = "complete.obs")
      
      # Energy
      obj_O3 <- weights$O3 * sum(abs(pcm - scm))
      
      # O1
      # Compute the break points (continuous sample quantiles)
      probs <- seq(0, 1, length.out = n_pts + 1)
      breaks <- lapply(covars[, id_num], stats::quantile, probs, na.rm = TRUE)
      
      # Count the number of points per marginal sampling strata and compare 
      # with the expected count
      sm_count <- lapply(1:length(id_num), function (i) 
        graphics::hist(sm[, id_num], breaks, plot = FALSE)$counts - 1)
      
      # Energy
      obj_O1 <- weights$O1 * sum(abs(sm_count))
    }
    
    if (any(covars_type == c("factor", "both")) {
      
      # O2
      
      # Compute the proportion of population points per marginal factor level
      pop_prop <- lapply(covars[, id_fac], function(x) table(x) / n_candi)
      
      # Compute the sample proportions
      sm_prop <- lapply(sm[, id_fac], function(x) table(x) / n_pts)
      
      # Compare the sample and population proportions
      sm_prop <- sapply(1:length(id_fac), function (i)
        sum(abs(sm_prop[[i]] - pop_prop[[i]])))
      
      # Energy (O2)
      obj_O2 <- weights$O2 * sum(sm_prop)
    }
    
    # Initial energy state
    if (covars_type == "both") {
      obj <- obj_O1 + obj_O2 + obj_O3
      energy0 <- data.frame(obj = obj, O1 = obj_O1, O2 = obj_O2, O3 = obj_O3)
      best_energy <- data.frame(obj = Inf, O1 = Inf, O2 = Inf, O3 = Inf)
    }
    if (covars_type == "numeric") {
      energy0 <- data.frame(obj = obj_O1 + obj_O3, O1 = obj_O1, O3 = obj_O3)
      best_energy <- data.frame(obj = Inf, O1 = Inf, O3 = Inf)
    }
    if (covars_type == "factor") {
      energy0 <- data.frame(obj = obj_O2)
      best_energy <- data.frame(obj = Inf)
    }

    # Other settings for the simulated annealing algorithm
    # ASR: How can we update these objects?
    # if (any(covars_type == c("numeric", "both")) {
      # old_scm <- scm
      # new_scm <- scm
      # best_scm <- scm
      # old_sm_count <- sm_count
      # new_sm_count <- sm_count
      # best_sm_count <- sm_count
    # }
    # if (any(covars_type == c("factor", "both")) {
      # old_sm_prop <- sm_prop
      # new_sm_prop <- sm_prop
      # best_sm_prop <- sm_prop
    # }
    old_sm <- sm
    new_sm <- sm
    best_sm <- sm
    count <- 0
    old_energy <- energy0
    if (progress) { 
      pb <- utils::txtProgressBar(min = 1, max = iterations, style = 3) 
    }
    time0 <- proc.time()

    # Begin the main loop
    for (k in 1:iterations) {
      
      # Plotting and jittering
      eval(.plot_and_jitter())
      
      # Update sample matrix
      # Recompute sample correlation matrix, marginal distribution, and
      # class proportions.
      # Recompute energy state
      new_sm[wp, ] <- covars[new_conf[wp, 1], ]
      
      if (any(covars_type == c("numeric", "both")) {
        
        # O3
        new_scm <- stats::cor(x = new_sm[, id_num], use = "complete.obs")
        obj_O3 <- weights$O3 * sum(abs(pcm - new_scm))
        
        # O1
        new_sm_count <- lapply(1:length(id_num), function (i)
          graphics::hist(new_sm[, id_num], breaks, plot = FALSE)$counts - 1)
        obj_O1 <- weights$O1 * sum(abs(new_sm_count))
        
      }
      if (any(covars_type == c("factor", "both")) {
        new_sm_prop <- lapply(sm[, id_fac], function(x) table(x) / n_pts)
        new_sm_prop <- sapply(1:length(id_fac), function (i)
          sum(abs(new_sm_prop[[i]] - pop_prop[[i]])))
        
        # Energy (O2)
        obj_O2 <- weights$O2 * sum(new_sm_prop)
      }
      
      if (covars_type == "both") {
        obj <- obj_O1 + obj_O2 + obj_O3
        new_energy <- 
          data.frame(obj = obj, O1 = obj_O1, O2 = obj_O2, O3 = obj_O3)
      }
      if (covars_type == "numeric") {
        new_energy <- 
          data.frame(obj = obj_O1 + obj_O3, O1 = obj_O1, O3 = obj_O3)
      }
      if (covars_type == "factor") {
        new_energy <- data.frame(obj = obj_O2)
      }
      
      # Evaluate the new system configuration
      random_prob <- ifelse(greedy, 1, stats::runif(1))
      actual_prob <- acceptance[[1]] * exp(-k / acceptance[[2]])
      if (track) { accept_probs[k] <- actual_prob }
      if (new_energy[1] <= old_energy[1]) {
        old_conf <- new_conf
        old_energy <- new_energy
        count <- 0
        old_sm <- new_sm
        # ASR: How can we update these objects?
        # if (any(covars_type == c("numeric", "both")) {
          # old_scm <- new_scm
          # old_sm_count <- new_sm_count
        # }
        # if (any(covars_type == c("factor", "both")) {
          # old_sm_prop <- new_sm_prop
        # }
      } else {
        if (new_energy[1] > old_energy[1] & random_prob <= actual_prob) {
          old_conf <- new_conf
          old_energy <- new_energy
          count <- count + 1
          old_sm <- new_sm
          # ASR: How can we update these objects?
          # if (any(covars_type == c("numeric", "both")) {
            # old_scm <- new_scm
            # old_sm_count <- new_sm_count
          # }
          # if (any(covars_type == c("factor", "both")) {
            # old_sm_prop <- new_sm_prop
          # }
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... p = ",
                random_prob, "\n")
          }
        } else {
          new_energy <- old_energy
          new_conf <- old_conf
          count <- count + 1
          new_sm <- old_sm
          # ASR: How can we update these objects?
          # if (any(covars_type == c("numeric", "both")) {
            # new_scm <- old_scm
            # new_sm_count <- old_sm_count
          # }
          # if (any(covars_type == c("factor", "both")) {
            # new_sm_prop <- old_sm_prop
          # }
          if (verbose) {
            cat("\n", count, "iteration(s) with no improvement... stops at",
                stopping[[1]], "\n")
          }
        }
      }
      
      # Best energy state
      if (track) energies[k, ] <- new_energy
      if (new_energy[1] < best_energy[1] / 1.0000001) {
        best_k <- k
        best_conf <- new_conf
        best_energy <- new_energy
        best_old_energy <- old_energy
        old_conf <- old_conf
        best_sm <- new_sm
        best_old_sm <- old_sm
        # ASR: How can we update these objects?
        # if (any(covars_type == c("numeric", "both")) {
          # best_scm <- new_scm
          # best_old_scm <- old_scm
          # best_sm_count <- new_sm_count
          # best_old_sm_count <- old_sm_count
        # }
        # if (any(covars_type == c("factor", "both")) {
          # best_sm_prop <- new_sm_prop
          # best_old_sm_prop <- old_sm_prop
        # }
      }
      
      # Freezing parameters
      if (count == stopping[[1]]) {
        if (new_energy[1] > best_energy[1] * 1.000001) {
          old_conf <- old_conf
          new_conf <- best_conf
          old_energy <- best_old_energy
          new_energy <- best_energy
          count <- 0
          new_sm <- best_sm
          old_sm <- best_old_sm
          # ASR: How can we update these objects?
          # if (any(covars_type == c("numeric", "both")) {
            # new_scm <- best_scm
            # old_scm <- best_old_scm
            # new_sm_count <- best_sm_count
            # old_sm_count <- best_old_sm_count
          # }
          # if (any(covars_type == c("factor", "both")) {
            # new_sm_prop <- best_sm_prop
            # old_sm_prop <- best_old_sm_prop
          # }
          cat("\n", "reached maximum count with suboptimal configuration\n")
          cat("\n", "restarting with previously best configuration\n")
          cat("\n", count, "iteration(s) with no improvement... stops at",
              stopping[[1]], "\n")
        } else {
          break
        }
      }
      if (progress) utils::setTxtProgressBar(pb, k)
    }

    # Prepare output
    eval(.prepare_output())
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS ##########################################
# candi: candidate locations
# covars: covariates
# use.coords: should the coordinates be used
.optimCLHScheck <-
  function (candi, covars, use.coords) {
    
    # covars
    if (is.vector(covars)) {
      if (use.coords == FALSE) {
        res <- "'covars' must have two or more columns"
        return (res)
      }
      if (nrow(candi) != length(covars)) {
        res <- "'candi' and 'covars' must have the same number of rows"
        return (res)
      }
    } else {
      if (nrow(candi) != nrow(covars)) {
        res <- "'candi' and 'covars' must have the same number of rows"
        return (res)
      }
    }
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE ############################
# This function is used to calculate the criterion value of CLHS.
# Aggregation is done using the weighted sum method.
.objCLHS <-
  function (sm, n.cov, weights, n.pts, pcm, scm, covars.type, pop.prop) {
    
    if (any(covars_type == c("numeric", "both")) {
      # O1
      obj_O1 <-
        .objO1(sm = sm, n.pts = n.pts, n.cov = n.cov, covars.type = covars.type)
      obj_O1 <- obj_O1 * weights$O1
      
      # O3
      obj_O3 <- .objO3(scm = scm, pcm = pcm)
      obj_O3 <- obj_O3 * weights$O3
    }
    
    if (any(covars_type == c("factor", "both")) {
    # O2
    obj_O2 <- 
      .objO2(sm = sm, n.pts = n.pts, n.cov = n.cov, pop.prop = pop.prop, 
             covars.type = covars.type)
    obj_O2 <- obj_O2 * weights$O2
    }
    
    # Prepare output
    if (covars_type == "both") {
      obj <- obj_O1 + obj_O2 + obj_O3
      res <- data.frame(obj = obj, O1 = obj_O1, O2 = obj_O2, O3 = obj_O3)
    }
    if (covars_type == "numeric") {
      res <- data.frame(obj = obj_O1 + obj_O3, O1 = obj_O1, O3 = obj_O3)
    }
    if (covars_type == "factor") {
      res <- data.frame(obj = obj_O2)
    }
    return (res)
  }
# CALCULATE OBJECTIVE FUNCTION VALUE ###########################################
#' @rdname optimCLHS
#' @export
objCLHS <-
  function (points, candi, covars, use.coords = FALSE, 
            weights = list(O1 = 1/3, O2 = 1/3, O3 = 1/3)) {
    
    # Check arguments
    check <- 
      .optimCLHScheck(candi = candi, covars = covars, use.coords = use.coords)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare 'covars' and create the starting sample matrix 'sm'
    eval(.prepare_acdc_covars())
    
    # Compute base data
    pcm <- .corCORR(obj = covars, covars.type = covars.type)
    scm <- .corCORR(obj = sm, covars.type = covars.type)
    pop_prop <- .strataACDC(n.pts = n_pts, strata.type = strata.type,
                            covars = covars, covars.type = covars.type)
    nadir <- .nadirACDC(n.pts = n_pts, n.cov = n_cov, n.candi = n_candi, 
                        pcm = pcm, nadir = nadir, candi = candi, 
                        covars = covars, pop.prop = pop_prop, 
                        covars.type = covars.type)
    utopia <- .utopiaACDC(utopia = utopia)
    
    # Compute the energy state
    energy <- .objACDC(sm = sm, pop.prop = pop_prop, covars.type = covars.type, 
                       weights = weights, pcm = pcm, scm = scm, n.pts = n_pts, 
                       n.cov = n_cov, utopia = utopia, nadir = nadir)
    return (energy)
  }
