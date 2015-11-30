# Version 1.0-2.9009 (2015-11-30)
* FIX: the computation of the number of point-pairs per lag-distance class in
  `optimPPL` was incorrect because it neglected the fact that, in a full 
  distance matrix, two points *a* and *b* form two pairs, i.e. *ab* and *ba*. 
  The mistake is due to the fact that we use `SpatialTools::dist1` to compute 
  the distance matrix instead of `stats::dist`.
* FEATURE: using a faster code to compute the number of points and point-pairs
  per lag-distance class in `optimPPL`.

# Version 1.0-2.9008 (2015-11-19)
* Fixed minor bugs.

# Version 1.0-2.9007 (2015-11-18)
* Improved the warning message printed when converting numeric covariates into
  factor covariates.
* Created a new `autofun` to check the number of accepted jitters in the first
  chain. If the number of accepted jitters is superior to the value passed to
  `schedule$initial.acceptance`, the process continues and a message is printed
  informing the proportion of jitters that have been accepted.
* Included scaling factors in two of the objective functions of `optimCLHS()`
  following the original Fortran code of Budiman Minasny.
* Use grey colours in plot with energy states; using only different line types
  was not enough to see the different lines -- using different colours makes it
  easier to see the differences among lines that represent different objective
  functions.
* Fixed minor bugs.

# Version 1.0-2.9006 (2015-11-17)
* The user can now chose the type of progress bar that should be used, with 
  options `"txt"`, for a text progress bar in the R console, `"tk"`, to put
  up a Tk progress bar widget, and `NULL` to omit the progress bar. A Tk 
  progress bar widget is useful when running ***spsann*** in parallel
  processors. The ***tcltk***-package is now a suggested package.
* Now we use grey colours to in the plot with the energy states.

# Version 1.0-2.9005 (2015-11-17)
* Solved NOTEs produced during CRAN check due to the use of functions from 
  default packages other than `base`, and due to examples that take more than 5
  seconds to run.

# Version 1.0-2.9004 (2015-11-16)
* Created a function to plot the optimized sample configuration (`plotOSC()`),
  with options to display the evolution of the energy state and/or the
  optimized sample configuration.
* The function used to compute the Pareto maximum and minimum values
  (`minmaxPareto()`) was optimized to be used with both ACDC and SPAN.

# Version 1.0-2.9003 (2015-11-15)
* Create a class (`OptimizedSampleConfiguration`) to store the output of
  `optim` functions.
  
# Version 1.0-2.9002 (2015-11-14)
* The trick included in the `optimMKV()`-function to avoid errors due to
  the LDLfactor error of the ***gstat***-package had to be reformulated. We are
  now using `try()` with a default value which is returned in case of error.

# Version 1.0-2.9001 (2015-11-13)
* A completelly new annealing schedule was implemented. The reason for this
  modification is that the former annealing schedule, which was based on the
  ***intamapInteractive***-package, showed to be ineficient during
  our tests. The new annealing schedule is the very simple and most-used
  schedule proposed by Kirkpatrick et al. (1983). We have also replaced the 
  acceptance criterion used in the ***intamapInteractive***-package with the
  well-known Metropolis criterion. This new implementation showed to be more
  efficient in our tests than our early implementation.
* Implementing a new annealing schedule and a new acceptance criterion required
  a moderate modification of the source code. Despite our efforts, it was not 
  possible to guarantee the compatibility with previous versions.
* A new function was created to set up the annealing schedule: 
  `scheduleSPSANN()`.
* We are now using a more elegant solution to jitter the sample points. It 
  consists of using a finite set of candidate locations that are seen by the
  algorithm as the centre of grid cells. In the first stage, we select a grid
  cell with replacement. In the second stage, we select a location within that
  grid cell using simple random sampling. This is the sampling method 
  implemented in the ***spcosa***-package.
* The documentation of all functions has been fine tuned.
* A trick was included in the `optimMKV()`-function to avoid errors due to
  the LDLfactor error of the ***gstat***-package.
* There also is a new function to compute the Pareto maximum and minimum values
  of the objective functions that compose a multi-objective optimization 
  problem (MOOP): `minmaxPareto()`.

# Version 1.0-2.9000 (2015-10-27)
* Now `x.max` and `y.max` are, by default, set to half of the maximum distance
  in the x- and y-coordinates of `candi`, respectively. In the same manner, the
  argument `cutoff` of `optimPPL()` is set, by default, to half of the diagonal
  of the rectangle of sides `x.max` and `y.max`.

# Version 1.0-2 (2015-10-25)
* Corrected a bug in `optimCORR()` that was causing the following error: Error
  in if (new_energy <= old_energy) { : missing value where TRUE/FALSE needed. 
  This bug used to affect `optimACDC()` and `optimSPAN()`.
* Created a version of the method proposed by Minasny and McBratney (2006),
  known as the conditioned Latin hypercube sampling (`optimCLHS`).
* Improved and updated the documentation of several functions.
* Improved the plotting functionality: each plot (the evolution of the
  energy state and the new system configuration) is now displayed in a separate
  device. This allows for a better visualization and allows the user to focus on
  a single plot if so s/he wishes.

# Version 1.0-1 (2015-07-30)
* Improved and updated documentation.
* ***gstat*** is not a dependence any more.
* Fixed breaks due to changes in dependencies (***pedometrics***).

# Version 1.0-0 (2015-07-14)
* Submission to CRAN.

# Version 0.0.0.9012 (2015-07-14)
* An auxiliary function (`objSPSANN()`) was created to retrieve the energy state
  of an optimized sample configuration (OSC) at a given point of the 
  optimization.
* Long examples are not run any more to avoid overload of `R CMD check`.
* The authors' list was updated with the respective roles.

# Version 0.0.0.9011 (2015-07-13)
* The documentation of all functions was significantly improved.
* Functions from default packages other than ***base*** are now imported to
  comply with the new change to the CRAN policy described at
  http://developer.r-project.org/blosxom.cgi/R-devel/NEWS/2015/06/29#n2015-06-29.

# Version 0.0.0.9010 (2015-07-05)
* Using `utils::globalVariables` to avoid the `R CMD check` note 
  `no visible binding for global variable [variable name]`. Source of the 
  solution: http://stackoverflow.com/a/12429344/3365410.
* The package ***fields*** is not a dependency any more.
* New default values were attributed to the following arguments: `plotit`, 
  `track`, `verbose`, and `iteration`. The first three were set to `FALSE`, 
  while the last was set to `100`.
* `optimSPAN()` and `objSPAN()` are now full operational.

# Version 0.0.0.9009 (2015-06-30)
* Several internal function were renamed using a pattern that includes the name
  of the respective objective function. For example, `.optimPPLcheck()` was 
  renamed as `.checkPPL()`, and `.getLagBreaks()` was renamed as `.lagsPPL()`.
  Note that the first part of the function name indicates what it does, while 
  the second indicates the objective function to which it applies. This 
  standardization is important to ease the construction of multi-objective
  optimization problems.

# Version 0.0.0.9008 (2015-06-29)
* Improvements in the family of ACDC, CORR, and DIST functions.
* Several pairs of internal function that were originally designed to deal with
  different types of covariates (factor and numeric) were merged. Now a single
  function does the job by using the key argument `covars.type`.
* New internal functions now enable building multi-objective optimization
  problems more easily. They have also allowed to clean-up/simplify the source 
  code.
* A new `autofun` was created to set-up the covariates (`covar`).

# Version 0.0.0.9007 (2015-06-12)
* The ***rgeos*** and ***plyr*** packages are not dependencies any more.
* The `boundary` of the spatial domain can now be estimated internally. The user 
  should use the ***rgeos*** package if a more precise `boundary` is needed.
* Now using a directory called 'R-autoFunction', where R code chunks that are 
  used in several functions of both families of `obj...()` and `optim...()`
  functions are included in individual files. These R code chunks are used to
  automatically build internal functions. Currently, R code chunks are
  used to check the arguments of the family of `optim...()` functions, prepare
  `points` and `candi`, set plotting options, estimate the `boundary`,
  prepare for jittering, plot and jitter, and prepare the output.
* BUG: the family of `obj...()` functions may not return the same criterion 
  value of the optimized sample configuration returned by the family of
  `optim...()` functions if the number of iterations used in the optimization 
  is equal to 100. The problem seems to disappear if a larger number of
  iterations is used.

# Version 0.0.0.9006 (2015-05-12)
* `spJitterFinite()` now tries to find an alternative point if the new point
  already is included in the sample. The number of tries is equal to the total
  number of points included in the sample. Because the more points we have, the
  more likely it is that the candidate point already is included in the sample.

# Version 0.0.0.9005 (2015-04-29)
* `spJitterFinite()` now returns the old point if the new point already is in
  the sample. This is to avoid an infinite loop at the end of the optimization
  when the objective function creates a cluster of points.

# Version 0.0.0.9004 (2015-04-24)
* New version of `optimACDC()`, including new argument definitions;
* In the multi-objective optimization problem case, now the graphical display 
  includes the many objective functions being optimized along with the utility
  function.

# Version 0.0.0.9003 (2015-04-20)
* Special version designed for the course on Spatial Sampling for Mapping, 
  22 - 24 April 2015, Wageningen University Campus, Wageningen, The Netherlands,
  Under the auspices of the Graduate School for Production Ecology and Resource 
  Conservation (PE&RC).

# Version 0.0.0.9002 (2015-04-19)
* new function to enable the user to define his/her own objective function;
* grammar check and enhanced documentation;

# Version 0.0.0.9001 (2015-02-24)
* new functions derived from `optimACDC()`: `optimDIST()` and `optimCORR()`;
* new objective function: mean/maximum kriging variance;
* review of the family of PPL functions;
* using function tailored argument checking.

# Version 0.0.0.9000 (2015-01-08)
* in-development package;
* importing functions from the package ***pedometrics***;
* preparing documentation.
