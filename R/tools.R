# Check spsann arguments #######################################################
.check_spsann_arguments <- 
  function (...) {
    expression(
      res <- NULL,
      aa <- c("points", "candi"),
      bb <- c(missing(points), missing(candi)), 
      if (any(bb)) {
        i <- which(bb == TRUE)
        res <- c("missing argument: '", aa[i], "'\n", sep = "")
      } else {
        # Argument 'candi'
        if (ncol(candi) != 2) {
          res <- c("'candi' must have two columns")
        } else {
          aa <- c("x", "y")
          bb <- colnames(candi)
          bb <- any(c(aa != bb) == TRUE)
          if (bb) {
            res <- c("'candi' must have two named columns: 'x' and 'y'")
          } else {
            # Argument 'iterations'
            if (!is.numint(iterations) || length(iterations) > 1) {
              res <- c("'iterations' must be an integer value")
            } else {
              # Argument 'acceptance'
              aa <- !is.list(acceptance)
              bb <- length(acceptance) != 2
              cc <- is.null(names(acceptance))
              dd <- !all(c(names(acceptance) == c("initial", "cooling")) == TRUE)
              if (aa || bb || cc || dd) {
                res <- paste("'acceptance' must be a list with two named ",
                             "sub-arguments: 'initial' and 'cooling'", sep = "")
              } else {
                # Argument 'stopping'
                aa <- !is.list(stopping)
                bb <- length(stopping) != 1
                cc <- is.null(names(stopping))
                dd <- !all(c(names(stopping) == "max.count") == TRUE)
                if (aa || bb || cc || dd) {
                  res <- paste("'stopping' must be a list with one named ",
                               "sub-argument: 'max.count'", sep = "")
                }
              }
            }
          }
        }
      }, 
      aa <- all(c(!missing(weights), !missing(utopia), !missing(nadir))), 
      if (aa) {
        aa <- !is.list(weights)
        bb <- is.null(names(weights))
        cc <- length(weights) < 2
        if (aa || bb || cc) {
          res <- c("'weights' must be a list with two or more named components")
        } else {
          aa <- sum(unlist(weights)) != 1
          if (aa) {
            res <- c("'weights' must sum to 1")
          } else {
            MOOP <- length(which(weights != 0))
            MOOP <- ifelse(MOOP > 1, TRUE, FALSE)
            COST <- ifelse(weights$COST == 0, FALSE, TRUE)
            
            # Argument 'utopia'
            if (MOOP) {
              aa <- !is.list(utopia)
              bb <- !length(utopia) == 1
              cc <- is.null(names(utopia))
              if (aa || bb || cc) {
                res <- c("'utopia' must be a list with one named component")
              } else {
                
                # Argument 'nadir'
                aa <- !is.list(nadir)
                if (aa) {
                  res <- c("'nadir' must be a list with named components")
                }
                aa <- names(nadir)
                if (length(aa) >= 2) {
                  if (length(aa) > 2) {
                    res <- c("you must choose a single 'nadir' option")
                  }
                } else {
                  if (aa == "sim" || aa == "seeds") { 
                    res <- c("the number of simulations and their seeds must be set")
                  }
                  if (aa == "abs") {
                    res <- c("sorry but 'nadir' cannot be calculated")
                  }
                }  
              }
            }
          }
        }
      },
      if (!is.null(res)) stop (res, call. = FALSE))
  }
# Set plotting options #########################################################
.plotting_options <-
  function (...) {
    expression(
      if (plotit) {
        par0 <- par()
        on.exit(suppressWarnings(par(par0)))
        if (missing(boundary)) {
          x <- by(candi[, 1], as.factor(candi[, 2]), 
                  FUN = range, simplify = FALSE)
          x <- do.call(rbind, x)
          d <- dist(x)
          d <- min(d[d > 0]) / 2
          x[, 1] <- x[, 1] - d
          x[, 2] <- x[, 2] + d
          y <- as.numeric(rownames(x))
          xy <- cbind(c(x[, 1], x[, 2]), rep(y, 2))
          
          y <- by(candi[, 2], as.factor(candi[, 1]), 
                  FUN = range, simplify = FALSE)
          y <- do.call(rbind, y)
          d <- dist(y)
          d <- min(d[d > 0]) / 2
          y[, 1] <- y[, 1] - d
          y[, 2] <- y[, 2] + d
          x <- as.numeric(rownames(y))
          yx <- cbind(rep(x, 2), c(y[, 1], y[, 2]))
          
          boundary <- SpatialPoints(unique(rbind(xy, yx)))
        }
      }
    )
  }
# Prepare for jittering ########################################################
.prepare_jittering <-
  function (...) {
    expression(
      if (missing(x.min) && missing(x.max) && missing(y.min) && missing(y.max)) {
        message("estimating 'x.min', 'x.max', 'y.min', and 'y.max' from 'candi'")
        x <- SpatialTools::dist1(as.matrix(candi[, "x"]))
        id <- x > 0
        x.min <- min(x[id])
        x.max <- max(x)
        y <- SpatialTools::dist1(as.matrix(candi[, "y"]))
        id <- y > 0
        y.min <- min(y[id])
        y.max <- max(y)
      },
      x_max0 <- x.max, y_max0 <- y.max
    )
  }
# Prepare points and candi #####################################################
.prepare_points <-
  function (...) {
    expression(
      if (!missing(candi)) {
        n_candi <- nrow(candi)
        candi <- as.matrix(cbind(id = 1:n_candi, candi))
      },
      if (is.integer(points) || is.numint(points)) {
        if (length(points) > 1) { # Integer vector
          points <- candi[points, ]
        }
        if (length(points) == 1) { # Integer value
          points <- sample(1:n_candi, points)
          points <- candi[points, ] 
        }
      } else { # Data frame of matrix
        points <- points
      },
      n_pts <- nrow(points), conf0 <- points, old_conf <- conf0
    )
  }
