% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minmaxPareto.R
\name{minmaxPareto}
\alias{minmaxPareto}
\title{Pareto minimum and maximum values}
\usage{
minmaxPareto(osc, candi, covars)
}
\arguments{
\item{osc}{A list of objects of class \code{OptimizedSampleConfiguration} (OSC). Each OSC of the list
must be named after the objective function with which it has been optimized. For example,
\code{osc = list(CORR = osc_corr, DIST = osc_dist)}.}

\item{candi}{Data frame (or matrix). The Cartesian x- and y-coordinates (in this order) of the
cell centres of a spatially exhaustive, rectangular grid covering the entire spatial sampling
domain. The spatial sampling domain can be contiguous or composed of disjoint areas and contain
holes and islands. \code{candi} provides the set of (finite) candidate locations inside the spatial
sampling domain for a point jittered during the optimization. Usually, \code{candi} will match the
geometry of the spatial grid containing the prediction locations, e.g. \code{newdata}
in \code{\link[gstat:krige]{gstat::krige()}}, \code{object} in \code{\link[raster:predict]{raster::predict()}}, and \code{locations} in \code{\link[geoR:krige.conv]{geoR::krige.conv()}}.}

\item{covars}{Data frame or matrix with the spatially exhaustive covariates in the columns.}
}
\value{
A data frame with the Pareto minimum and maximum values.
}
\description{
Compute the minimum and maximum attainable values of the objective functions that compose a
multi-objective combinatorial optimization problem.
}
\details{
\subsection{Multi-objective combinatorial optimization problems}{
A method of solving a multi-objective combinatorial optimization problem (MOCOP) is to aggregate
the objective functions into a single \emph{utility function}. In \strong{spsann}, the aggregation
is performed using the \emph{weighted sum method}, which incorporates in the weights the
preferences of the user regarding the relative importance of each objective function.

The weighted sum method is affected by the relative magnitude of the different function values.
The objective functions implemented in \strong{spsann} have different units and orders of magnitude.
The consequence is that the objective function with the largest values may have a numerical
dominance during the optimization. In other words, the weights may not express the true
preferences of the user, resulting that the meaning of the utility function becomes unclear
because the optimization will favour the objective function which is numerically dominant.

A reasonable solution to avoid the numerical dominance of any objective function is to scale the
objective functions so that they are constrained to the same approximate range of values. Several
function-transformation methods can be used for this end and \strong{spsann} has four of them
available.

The \emph{upper-lower-bound approach} requires the user to inform the maximum (nadir point) and
minimum (utopia point) absolute function values. The resulting function values will always range
between 0 and 1.

The \emph{upper-bound approach} requires the user to inform only the nadir point, while the utopia
point is set to zero. The upper-bound approach for transformation aims at equalizing only the
upper bounds of the objective functions. The resulting function values will always be smaller
than or equal to 1.

In most cases, the absolute maximum and minimum values of an objective function cannot be
calculated exactly. If the user is uncomfortable with guessing the nadir and utopia points, there
an option for using \emph{numerical simulations}. It consists of computing the function value for many
random system configurations. The mean function value obtained over multiple simulations is used
to set the nadir point, while the the utopia point is set to zero. This approach is similar to
the upper-bound approach, but the function values will have the same orders of magnitude only at
the starting point of the optimization. Function values larger than one are likely to occur
during the optimization. We recommend the user to avoid this approach whenever possible because
the effect of the starting configuration on the optimization as a whole usually is insignificant
or arbitrary.

The \emph{upper-lower-bound approach} with the minimum and maximum \emph{attainable} values of the
objective functions that compose the MOCOP, also known as the \emph{Pareto minimum and maximum
values}, is the most elegant solution to scale the objective functions. However, it is the most
time consuming. It works as follows:
\enumerate{
\item Optimize a sample configuration with respect to each objective function that composes the
MOCOP;
\item Compute the function value of every objective function that composes the MOCOP for every
optimized sample configuration;
\item Record the minimum and maximum absolute function values attained for each objective function
that composes the MOCOP -- these are the Pareto minimum and maximum.
}

For example, consider ACDC, a MOCOP composed of two objective functions: CORR and FREQ. The
minimum absolute attainable value of CORR is obtained when the sample configuration is optimized
with respect to only CORR, i.e. when the evaluator and generator objective functions are the same
(see the intersection between the second line and second column in the table below). This is
the Pareto minimum of CORR. It follows that the maximum absolute attainable value of CORR is
obtained when the sample configuration is optimized with regard to only FREQ, i.e. when the
evaluator function is difference from the generator function (see the intersection between the
first row and the second column in the table below). This is the Pareto maximum of CORR. The same
logic applies for finding the Pareto minimum and maximum of FREQ.

\tabular{rll}{
\emph{Evaluator} \tab \emph{Generator} \tab      \cr
\tab FREQ             \tab CORR \cr
FREQ             \tab 0.5              \tab 8.6  \cr
CORR             \tab 6.4              \tab 0.3  \cr
}

}
}
\examples{
#####################################################################
# NOTE: The settings below are unlikely to meet your needs.         #
#####################################################################
if (interactive() & require(sp)) {
  # This example takes more than 5 seconds to run
  data(meuse.grid, package = "sp")
  # General (greedy) cooling schedule
  schedule <- scheduleSPSANN(
    initial.acceptance = c(0.01, 0.99), chains = 100,
    x.max = 1540, y.max = 2060, x.min = 0,
    y.min = 0, cellsize = 40)
  # CORR
  schedule$initial.temperature <- 0.1
  set.seed(2001)
  osc_corr <- optimCORR(
    points = 10, candi = meuse.grid[, 1:2],
    covars = meuse.grid[, c(1, 2)],
    schedule = schedule)
  # FREQ
  schedule$initial.temperature <- 0.1
  set.seed(2001)
  osc_dist <- optimDIST(
    points = 10, candi = meuse.grid[, 1:2],
    covars = meuse.grid[, c(1, 2)],
    schedule = schedule)
  # PPL
  schedule$initial.temperature <- 0.1
  set.seed(2001)
  osc_ppl <- optimPPL(
    points = 10, candi = meuse.grid[, 1:2],
    schedule = schedule)
  # MSSD
  schedule$initial.temperature <- 0.1
  set.seed(2001)
  osc_mssd <- optimMSSD(
    points = 10, candi = meuse.grid[, 1:2],
    schedule = schedule)
  # Pareto
  pareto <- minmaxPareto(
    osc = list(
      DIST = osc_dist,
      CORR = osc_corr,
      PPL = osc_ppl,
      MSSD = osc_mssd),
    candi = meuse.grid[, 1:2],
    covars = meuse.grid[, c(1, 2)])
  round(pareto, 4)
}
}
\references{
Arora, J. \emph{Introduction to optimum design}. Waltham: Academic Press, p. 896, 2011.

Marler, R. T.; Arora, J. S. Survey of multi-objective optimization methods for engineering.
\emph{Structural and Multidisciplinary Optimization}, v. 26, p. 369-395, 2004.

Marler, R. T.; Arora, J. S. Function-transformation methods for multi-objective optimization.
\emph{Engineering Optimization}, v. 37, p. 551-570, 2005.

Marler, R. T.; Arora, J. S. The weighted sum method for multi-objective optimization: new
insights. \emph{Structural and Multidisciplinary Optimization}, v. 41, p. 853-862, 2009.
}
\seealso{
\code{\link[=optimACDC]{optimACDC()}}, \code{\link[=SPAN]{SPAN()}}
}
\author{
Alessandro Samuel-Rosa \email{alessandrosamuelrosa@gmail.com}
}
