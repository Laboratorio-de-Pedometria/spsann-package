% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimMKV.R
\name{optimMKV}
\alias{optimMKV}
\alias{objMKV}
\alias{MKV}
\title{Optimization of sample configurations for spatial interpolation (II)}
\usage{
optimMKV(
  points,
  candi,
  covars,
  eqn,
  vgm,
  krige.stat = "mean",
  ...,
  schedule,
  plotit = FALSE,
  track = FALSE,
  boundary,
  progress = "txt",
  verbose = FALSE
)

objMKV(points, candi, covars, eqn, vgm, krige.stat = "mean", ...)
}
\arguments{
\item{points}{Integer value, integer vector, data frame (or matrix), or list. The number of
sampling points (sample size) or the starting sample configuration. Four options are available:
\itemize{
\item Integer value. The required number of sampling points (sample size). The sample configuration
used to start the optimization will consist of grid cell centres of \code{candi} selected using
simple random sampling, i.e. \code{\link[base:sample]{base::sample()}} with \code{x = 1:nrow(candi)} and \code{size = points}.
\item Integer vector. A set of row indexes between one (1) and \code{nrow(candi)}. These row indexes
identify the grid cell centres of \code{candi} that will form the starting sample configuration for
the optimization. The length of the integer vector, \code{length(points)}, is the sample size.
\item Data frame (or matrix). The Cartesian x- and y-coordinates (in this order) of the starting
sample configuration.
\item List. An object with two named sub-arguments:
\itemize{
\item \code{fixed} An integer vector or data frame (or matrix) specifying an existing sample
configuration (see options above). This sample configuration is kept as-is (fixed) during
the optimization and is used only to compute the objective function values.
\item \code{free} An integer value, integer vector, data frame or matrix (see options above) specifying
the (number of) sampling points to add to the existing sample configuration. These new
sampling points are free to be moved around (jittered) during the optimization.
}
}

Most users will want to set an integer value simply specifying the required sample size. Using
an integer vector or data frame (or matrix) will generally be helpful to users willing to
evaluate starting sample configurations, test strategies to speed up the optimization, and
fine-tune or thin an existing sample configuration. Users interested in augmenting a possibly
existing real-world sample configuration or fine-tuning only a subset of the existing sampling
points will want to use a list.}

\item{candi}{Data frame (or matrix). The Cartesian x- and y-coordinates (in this order) of the
cell centres of a spatially exhaustive, rectangular grid covering the entire spatial sampling
domain. The spatial sampling domain can be contiguous or composed of disjoint areas and contain
holes and islands. \code{candi} provides the set of (finite) candidate locations inside the spatial
sampling domain for a point jittered during the optimization. Usually, \code{candi} will match the
geometry of the spatial grid containing the prediction locations, e.g. \code{newdata}
in \code{\link[gstat:krige]{gstat::krige()}}, \code{object} in \code{\link[raster:predict]{raster::predict()}}, and \code{locations} in \code{\link[geoR:krige.conv]{geoR::krige.conv()}}.}

\item{covars}{Data frame or matrix with the covariates in the columns. The number of rows of \code{covars} must
match exactly that of \code{candi} -- or \code{eval.grid}, in case a coarser evaluation grid is used.}

\item{eqn}{Formula string that defines the dependent variable \code{z} as a linear function of the independent
variables (covariates) contained in \code{covars}. See the argument \code{formula} in the function
\verb{\link[gstat]\{krige\}} for more information.}

\item{vgm}{Object of class \code{variogramModel}. See the argument \code{model} in the function \verb{\link[gstat]\{krige\}}
for more information.}

\item{krige.stat}{Character value defining the statistic that should be used to summarize the kriging
prediction error variance. Available options are \code{"mean"} and \code{"max"} for the mean and maximum kriging
prediction error variance, respectively. Defaults to \code{krige.stat = "mean"}.}

\item{...}{further arguments passed to \verb{\link[gstat]\{krige\}}. (Advanced users only!)}

\item{schedule}{List with named sub-arguments setting the control parameters of the annealing
schedule. See \code{\link[=scheduleSPSANN]{scheduleSPSANN()}}.}

\item{plotit}{(Optional) Logical for plotting the evolution of the optimization. Plot updates
occur at each ten (10) spatial jitters. Defaults to \code{plotit = FALSE}. The plot includes two
panels:
\enumerate{
\item The first panel depicts the changes in the objective function value (y-axis) with the
annealing schedule (x-axis). The objective function values should be high and variable at the
beginning of the optimization (panel's top left). As the optimization proceeds, the objective
function values should gradually transition to a monotone decreasing behaviour till they
become virtually constant. The objective function values constancy suggests the end of the
optimization (panel's bottom right).
\item The second panel shows the starting (grey circles) and current spatial sample configuration
(black dots). Black crosses indicate the fixed (existing) sampling points when a spatial
sample configuration is augmented. The plot shows the starting sample configuration to assess
the effects on the optimized spatial sample configuration: the latter generally should be
independent of the first. The second panel also shows the maximum possible spatial jitter
applied to a sampling point in the Cartesian x- (x-axis) and y-coordinates (y-axis).
}}

\item{track}{(Optional) Logical value. Should the evolution of the energy state be recorded and
returned along with the result? If \code{track = FALSE} (the default), only the starting and ending
energy states return along with the results.}

\item{boundary}{(Optional) An object of class SpatialPolygons (see sp::SpatialPolygons()) with
the outer and inner limits of the spatial sampling domain (see \code{candi}). These SpatialPolygons
help depict the spatial distribution of the (starting and current) sample configuration inside
the spatial sampling domain. The outer limits of \code{candi} serve as a rough \code{boundary} when
\code{plotit = TRUE}, but the SpatialPolygons are missing.}

\item{progress}{(Optional) Type of progress bar that should be used, with options \code{"txt"}, for a
text progress bar in the R console, \code{"tk"}, to put up a Tk progress bar widget, and \code{NULL} to
omit the progress bar. A Tk progress bar widget is useful when using parallel processors.
Defaults to \code{progress = "txt"}.}

\item{verbose}{(Optional) Logical for printing messages about the progress of the optimization.
Defaults to \code{verbose = FALSE}.}
}
\value{
\code{optimMKV} returns an object of class \code{OptimizedSampleConfiguration}: the optimized sample
configuration with details about the optimization.

\code{objMKV} returns a numeric value: the energy state of the sample configuration -- the objective
function value.
}
\description{
Optimize a sample configuration for spatial interpolation with a 'known' linear mixed model, e.g. universal
(external drift) kriging and regression-kriging with a linear regression model. A criterion is defined so
that the sample configuration minimizes the mean or maximum kriging prediction error variance (\bold{MKV}).
}
\details{
\subsection{Generating mechanism}{
There are multiple mechanism to generate a new sample configuration out of an existing one. The
main step consists of randomly perturbing the coordinates of a single sample, a process known as
\sQuote{jittering}. These mechanisms can be classified based on how the set of candidate
locations for the samples is defined. For example, one could use an \emph{infinite} set of candidate
locations, that is, any location in the spatial domain can be selected as a new sample location
after a sample is jittered. All that is needed is a polygon indicating the boundary of the
spatial domain. This method is more computationally demanding because every time an existing
sample is jittered, it is necessary to check if the new sample location falls in spatial domain.

Another approach consists of using a \emph{finite} set of candidate locations for the samples. A
finite set of candidate locations is created by discretising the spatial domain, that is,
creating a fine (regular) grid of points that serve as candidate locations for the jittered
sample. This is a less computationally demanding jittering method because, by definition, the
new sample location will always fall in the spatial domain.

Using a finite set of candidate locations has two important inconveniences. First, not all
locations in the spatial domain can be selected as the new location for a jittered sample.
Second, when a sample is jittered, it may be that the new location already is occupied by another
sample. If this happens, another location has to be iteratively sought for, say, as many times as
the size of the sample configuration. In general, the larger the size of the sample
configuration, the more likely it is that the new location already is occupied by another sample.
If a solution is not found in a reasonable time, the the sample selected to be jittered is kept
in its original location. Such a procedure clearly is suboptimal.

\strong{spsann} uses a more elegant method which is based on using a finite set of candidate locations
coupled with a form of \emph{two-stage random sampling} as implemented in \code{\link[spcosa:spsample-methods]{spcosa::spsample()}}.
Because the candidate locations are placed on a finite regular grid, they can be taken as the
centre nodes of a finite set of grid cells (or pixels of a raster image). In the first stage, one
of the \dQuote{grid cells} is selected with replacement, i.e. independently of already being
occupied by another sample. The new location for the sample chosen to be jittered is selected
within that \dQuote{grid cell} by simple random sampling. This method guarantees that virtually
any location in the spatial domain can be selected. It also discards the need to check if the new
location already is occupied by another sample, speeding up the computations when compared to the
first two approaches.
}

\subsection{Annealing schedule}{
The \emph{search graph} corresponds to the set of effective candidate locations for a sample location
selected to be jittered. The size of the search graph, i.e. area within which a sample location
can be moved around, is related to the concept of \emph{temperature}. A larger search graph is
equivalent to higher temperatures, which potentially result in more movement -- or
\sQuote{agitation} -- of the set of sample locations.

The current version of the \pkg{spsann}-package uses a linear cooling schedule which depends upon
the number of jitters to control the size of the search graph. The equations are

\eqn{x_max = x_max0 - (chains_i / chains) * (x_max0 - x_min) + x_cellsize + x_min0}

and

\eqn{y_max = y_max0 - (chains_i / chains) * (y_max0 - y_min) + y_cellsize + y_min0},

where $x_max0$ and $y_max0$ are the maximum allowed shifts in the x- and y-coordinates in the
first chain, $x_min$ and $y_min$ are the minimum required shifts in the x- and y-coordinates,
$x_max$ and $y_max$ are the maximum allowed shifts in the x- and y-coordinates during the next
chain, $chains$ and $chain_i$ are the total and current chains, and $x_cellsize$ and $y_cellsize$
are the grid spacing in the x- and y-coordinates. Because $x_cellsize$ and $y_cellsize$ can be
equal to zero when a finite set of candidate locations is used, $x_min0$ and $y_min0$ are the
maximum nearest neighbour distance in the x- and y-coordinates between candidate locations.
}
}
\note{
\subsection{Distance between two points}{
\strong{spsann} always computes the distance between two locations (points) as the
\href{https://en.wikipedia.org/wiki/Euclidean_distance}{Euclidean distance} between them. This
computation requires the optimization to operate in the two-dimensional Euclidean space, i.e. the
coordinates of the sample, candidate and evaluation locations must be Cartesian coordinates,
generally in metres or kilometres. \strong{spsann} has no mechanism to check if the coordinates are
Cartesian: you are the sole responsible for making sure that this requirement is attained.
}

This function is based on the method originally proposed by Heuvelink, Brus and de Gruijter (2006) and
implemented in the R-package \pkg{intamapInteractive} by Edzer Pebesma and Jon Skoien.
}
\examples{
#####################################################################
# NOTE: The settings below are unlikely to meet your needs.         #
#####################################################################
\dontrun{
data(meuse.grid, package = "sp")
candi <- meuse.grid[1:1000, 1:2]
covars <- as.data.frame(meuse.grid)[1:1000, ]
vgm <- gstat::vgm(psill = 10, model = "Exp", range = 500, nugget = 8)
schedule <- scheduleSPSANN(
  initial.temperature = 10, chains = 1, x.max = 1540, y.max = 2060, 
  x.min = 0,  y.min = 0, cellsize = 40)
set.seed(2001)
res <- optimMKV(
  points = 10, candi = candi, covars = covars, eqn = z ~ dist, 
  vgm = vgm, schedule = schedule)
data.frame(
  expected = 15.37137,
  objSPSANN = objSPSANN(res),
  objMKV = objMKV(
    points = res, candi = candi, covars = covars, eqn = z ~ dist, vgm = vgm)
)
}
}
\references{
Brus, D. J.; Heuvelink, G. B. M. Optimization of sample patterns for universal kriging of environmental
variables. \emph{Geoderma}. v. 138, p. 86-95, 2007.

Heuvelink, G. B. M.; Brus, D. J.; de Gruijter, J. J. Optimization of sample configurations for digital
mapping of soil properties with universal kriging. In: Lagacherie, P.; McBratney, A. & Voltz, M. (Eds.)
\emph{Digital soil mapping - an introductory perspective}. Elsevier, v. 31, p. 137-151, 2006.
}
\author{
Alessandro Samuel-Rosa \email{alessandrosamuelrosa@gmail.com}
}
\concept{simulated annealing}
\concept{spatial interpolation}
\keyword{iteration}
\keyword{optimize}
\keyword{spatial}
