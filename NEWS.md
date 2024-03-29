# NEWS

## spsann 2.2.3

* __Improvement.__ Package documentation is expanded, with more information for inexperienced users.

## spsann 2.2.2

* __Improvement.__ Package documentation is expanded, with more information for inexperienced users.

## spsann 2.2.1

* __Improvement.__ Package documentation is expanded, with more information for inexperienced users.

## spsann 2.2.0.9009

* __Improvement.__ PPL-related functions. Improved argument checking using `match.call()`, the
  code now being cleaner. `schedule` does not have a default value any more. This is in line with
  version 2.2.0.9002.

## spsann 2.2.0.9008

* __Improvement 1.__ `optimPPL()`: the optimized sample configuration now exports both target
  (`target.distribution`) and final (`final.distribution`) distribution of points or point-pairs
  per lag-distance class.
* __Improvement 2.__ `objPPL()`: for a PPL-optimized sample configuration, now returns the stored
  objective function value - `objSPSANN()` may become useless.
* __Improvement 3.__ `countPPL()`: for a PPL-optimized sample configuration, now returns the stored
  final distribution of points or point-pairs per lag-distance class.

## spsann 2.2.0.9007

* __Bug fix.__ `minmaxPareto()`: `'cellsize' must be a numeric value greater than or equal to zero`.
* __Improvement.__ `minmaxPareto()`: updates example in documentation.

## spsann 2.2.0.9006 (2020-08-06)

* Improvements:
  * Now using a qualitative, colour-blind safe colour palette based on ColorBrewer (<https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=4>)

## spsann 2.2.0.9005 (2020-07-24)

* Improvements:
  * Function `scheduleSPSANN`: Argument `initial.acceptance` now requires
    a numeric vector of length two defining the minimum and maximum initial
    acceptance probability, with default `c(0.95, 0.99)`. This was needed
    so that users have more confidence choosing the `initial.temperature`.
    Some were using unnecessarily high initial temperature values.
  * Improves documentation to accommodate changes.

## spsann 2.2.0.9004 (2019-10-01)

* Improvements:
  * Cleans up code used to estimate jittering parameters.
  * Improves documentation to accommodate latest changes/fixes.
  * Argument `stopping` of function `scheduleSPSANN()` is now set to be equal to fifty percent the
    maximum number of chains by using `ceiling(chains * 0.5)`.

## spsann 2.2.0.9003 (2019-09-30)

* Bug fix:
  * Avoid points being stuck in the same location at the end of the optimization when using a
    finite set of candidate locations. This is done by computing 'x_min0' and 'y_min0'. These are
    the maximum distance to the nearest neighbouring candidate location in the x- and y-coordinates.

## spsann 2.2.0.9002 (2019-09-15)

* Improvements:
  * Testing `optimMKV()` as a greedy algorithm with (1) a too small neighbourhood size and (2) with
    a neighbourhood set using `nmax`. The first throws an error -- `skipped 'singular matrix' error
    in 'krige'-function` -- due to the too small number of samples in the neighbourhood.
  * Argument `eval.grid` temporarily removed from functions `optimMKV()` and `objMKV()` and flagged
    as experimental in functions `optimMSSD()` and `objMSSD()`. The reason is that the
    implementation for `optimMKV()` and `objMKV()` is less trivial because of the presence of
    covariates that need to match `eval.grid`.
  * Argument `cellsize` of function `scheduleSPSANN` is not estimated from `candi` any more. This
    is because `candi` can also be an existing irregular sample configuration and inexperienced
    users can be unaware of the need to set `cellsize = 0` in this case. Only `x.max` and `y.max`
    are still estimated from `candi`. The user now has to inform `cellsize` manually and this
    should speed up computations.

## spsann 2.2.0.9001 (2019-09-14)

* Improvements:
  * Functions `optimMSSD()`, `objMSSD()`, `optimMKV()` and `objMKV()` have a new argument called `eval.grid`.
    With this new argument, one can use a grid that is coarser than `candi` to compute the objective function 
    value (energy state). This is important to speed up computations in large areas. With this new argument one
    can also thin an existing sample configuration. To do this, the existing sample configuration has to be 
    passed to `candi`. If one passes a matrix to `eval.grid`, then jittering parameters are estimated from 
    `eval.grid` instead of `candi`.
  * A note has been added to the examples of all functions. Is says: 'The settings below are unlikely to meet
    your needs.' This is an attempt to warn inexperienced users of the need to study each function and set
    function arguments according to their use cases.
  * Argument `eqn` of function `optimMKV()` and `objMKV()` has no default any more. The previous default --
    `eqn = z ~ 1` -- represented an ordinary kriging system, while the function was primarily designed for 
    universal kriging. The default was useful as an example, but has been causing confusion among inexperienced
     users, thus the decision for its removal.
  * Argument `stopping` of function `scheduleSPSANN()` is now set to be equal to ten percent the maximum
    number of chains by using `ceiling(chains * 0.1)`. The previous default value of `stopping = 10` was being
    misused by inexperienced users, thus the decision for its modification. With the new default settings,
    `chain = 500` and thus `stopping = 50`, which should give more stable results.
  * The documentation of various functions has been revised and improved wherever necessary.

## spsann 2.2.0.9000 (2019-05-06)

* Improvement:
    - Adds argument `verbose` to `spJitter()`. The user decides if s/he wants to be warned when there is no
      candidate location left in the neighbourhood and the sample point is kept in its original location.

## spsann 2.2.0 (2019-04-28)

The new version of the __spsann__ package includes some bug fixes and a few modifications. Users now can 
choose how `optimCLHS` computes objective function values: as in the original paper or as in the FORTRAN 
implementation. Users now also must inform the `weights` passed to `optimCLHS` as to guarantee that s/he is 
aware of what s/he is doing. The same apples to other functions that deal with multi-objective optimization 
problems: `optimACDC` and `optimSPAM`. Another important modification in the current version of __spsann__ is 
the possibility to use a finite set of candidate locations by setting `cellsize = 0`. This is useful when 
optimizing sample points only in the feature space and should reduce the computation time needed to find the
solution.

## spsann 2.1.0.9005 (2019-04-26)

* Accounts for the fact that the __[clhs](https://CRAN.R-project.org/package=clhs)__ package 
  by Pierre Roudier is not on CRAN any more.
* Corrects a bug in `plot.OptimizedSampleConfiguration` related to the selection of the information to be 
  displayed.

## spsann 2.1.0.9004 (2019-04-10)

* Argument `weights` of `optimCLHS`, `optimACDC` and `optimSPAM` now is mandatory. The user is now required to 
  set the weights as to guarantee that s/he is aware of what s/he is doing.
* Corrects bug related to the naming conventions used with `data.frame`s that store objective function
  values.
* Allow using a finite set of candidate locations for jittering sample points by setting `cellsize = 0`. When
  this is done, __spsann__ now checks for neighbouring candidate locations already included in the sample as to
  avoid duplicated sampling points.

## spsann 2.1.0.9003 (2018-07-19)

* Implements modifications (alternatives) to the way how `optimCLHS` computes objective function values.

## spsann 2.1.0.9002 (2018-07-18)

* Implements modifications (alternatives) to the way how `optimCLHS` computes objective function values.

## spsann 2.1.0.9001 (2018-07-17)

* Updates package documentation.

## spsann 2.1.0.9000 (2018-07-16)

* Adds new badges to README.md: package version and project status.
* Adds details and improves DESCRIPTION file; authors are reordered based on contributions.
* Implements modifications (alternatives) to the way how `optimCLHS` computes objective function values.

## spsann 2.1-0 (2017-06-23)

Now __spsann__ can be used to augment an existing sample configuration, that is, add new sampling points
to a spatial sample configuration generated using __spsann__ or any other means. To do so, when using one 
of the functions from the family of `optim...()` functions, the user must pass to the function argument 
`points` an object of class `list` containing two named sub-arguments: `fixed`, a matrix with the 
coordinates of the existing sample configuration -- kept fixed during the optimization --, and `free`,
the number of sample points that should be added to the existing sample configuration -- free to move around
during the optimization.

## spsann 2.0-0 (2016-03-14)

This is a major release of package __spsann__ that includes several conceptual changes. Despite our efforts, 
it was not possible to guarantee the compatibility with previous versions. We have decided not to deprecate 
functions and function arguments because (1) this would require deprecating a lot of code and (2) you should 
first read the updated package documentation to understand the conceptual changes that we have made before you
start using it. This is a summary of the changes:
* A completely new annealing schedule was implemented. The reason for this modification is that the former
  annealing schedule showed to be inefficient during our tests. The new annealing schedule is the very simple 
  and most-used schedule proposed by Kirkpatrick et al. (1983). We have also replaced the acceptance criterion 
  with the well-known Metropolis criterion. This new implementation showed to be more efficient in our tests 
  than our early implementation. Setting up this new annealing schedule is done using the new function 
  `scheduleSPSANN`.
* A more elegant solution to jitter the sample points was implemented. It consists of using a finite set of
  candidate locations that are seen by the algorithm as the centre of grid cells. In the first stage, we 
  select a grid cell with replacement. In the second stage, we select a location within that grid cell using
  simple random sampling. This guarantees that any location in the sampling region is a candidate location 
  for the jittered sample point.
* Solving multi-objective combinatorial optimization problems (MOCOP) has become easier with the creation of 
  the new function `minmaxPareto()`. This function computes the Pareto maximum and minimum values of the 
  objective functions that compose the MOCOP needed to scale the objective functions to the same approximate
  range of values.
* The user can now chose to follow the progress of the optimization using a text progress bar in the R console
  or a Tk progress bar widget. A Tk progress bar widget is useful when running __spsann__ in parallel
  processors.
* The output of the optimization is now stored in an object of class `OptimizedSampleConfiguration`. This 
  object contains three slots. The first (`points`) holds the coordinates of the optimized sample 
  configuration. The second, `spsann`, stores information about the settings used with the spatial simulated
  annealing algorithm. The third, `objective`, holds the settings used with the chosen objective function. 
  Methods were implemented to retrieve information from the new class, as well as producing plots of the 
  optimized sample configuration.
* Package documentation was expanded and adapted to cope with the conceptual changes that were made. It also 
  includes a vignette that gives a short description of the package and its structure, as well as presents 
  a few examples on how to use the package. It is strongly recommended to read the new package documentation 
  and the accompanying vignette before you start using the package. 
* Finally, bugs were fixed, warning messages were improved, and a faster code was implemented whenever 
  possible.

## spsann 1.0-2.9013 (2016-03-14)

* `devel` branch was merged into `master` branch.

## spsann 1.0-2.9012 (2016-03-14)

* Package documentation was expanded. It now includes a vignette that gives a short description of
  the package and its structure. The vignette also contains a few examples on how to use the
  package. `knitr` is the engine used to produce the package vignette.

## spsann 1.0-2.9011 (2016-03-11)

* Documentation was expanded.

## spsann 1.0-2.9010 (2016-03-10)

* Created S3 methods for extracting the objective function value and plotting an object of class
  `OptimizedSampleConfiguration`.
* The class `OptimizedSampleConfiguration` is no longer exported.
* Documentation was expanded.

## spsann 1.0-2.9009 (2015-11-30)

* FIX: the computation of the number of point-pairs per lag-distance class in
  `optimPPL` was incorrect because it neglected the fact that, in a full 
  distance matrix, two points *a* and *b* form two pairs, i.e. *ab* and *ba*. 
  The mistake is due to the fact that we use `SpatialTools::dist1` to compute 
  the distance matrix instead of `stats::dist`.
* FEATURE: using a faster code to compute the number of points and point-pairs
  per lag-distance class in `optimPPL`.

## spsann 1.0-2.9008 (2015-11-19)

* Fixed minor bugs.

## spsann 1.0-2.9007 (2015-11-18)

* Improved the warning message printed when converting numeric covariates into factor covariates.
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

## spsann 1.0-2.9006 (2015-11-17)

* The user can now chose the type of progress bar that should be used, with 
  options `"txt"`, for a text progress bar in the R console, `"tk"`, to put
  up a Tk progress bar widget, and `NULL` to omit the progress bar. A Tk 
  progress bar widget is useful when running __spsann__ in parallel
  processors. The __tcltk__-package is now a suggested package.
* Now we use grey colours to in the plot with the energy states.

## spsann 1.0-2.9005 (2015-11-17)

* Solved NOTEs produced during CRAN check due to the use of functions from 
  default packages other than `base`, and due to examples that take more than 5
  seconds to run.

## spsann 1.0-2.9004 (2015-11-16)

* Created a function to plot the optimized sample configuration (`plotOSC()`), with options to
  display the evolution of the energy state and/or the optimized sample configuration.
* The function used to compute the Pareto maximum and minimum values (`minmaxPareto()`) was
  optimized to be used with both ACDC and SPAN.

## spsann 1.0-2.9003 (2015-11-15)

* Create a class (`OptimizedSampleConfiguration`) to store the output of `optim` functions.

## spsann 1.0-2.9002 (2015-11-14)

* The trick included in the `optimMKV()`-function to avoid errors due to
  the LDLfactor error of the __gstat__-package had to be reformulated. We are
  now using `try()` with a default value which is returned in case of error.

## spsann 1.0-2.9001 (2015-11-13)

* A completely new annealing schedule was implemented. The reason for this
  modification is that the former annealing schedule, which was based on the
  __intamapInteractive__-package, showed to be inefficient during
  our tests. The new annealing schedule is the very simple and most-used
  schedule proposed by Kirkpatrick et al. (1983). We have also replaced the 
  acceptance criterion used in the __intamapInteractive__-package with the
  well-known Metropolis criterion. This new implementation showed to be more
  efficient in our tests than our early implementation.
* Implementing a new annealing schedule and a new acceptance criterion required
  a moderate modification of the source code. Despite our efforts, it was not 
  possible to guarantee the compatibility with previous versions.
* A new function was created to set up the annealing schedule: `scheduleSPSANN()`.
* We are now using a more elegant solution to jitter the sample points. It 
  consists of using a finite set of candidate locations that are seen by the
  algorithm as the centre of grid cells. In the first stage, we select a grid
  cell with replacement. In the second stage, we select a location within that
  grid cell using simple random sampling. This is the sampling method 
  implemented in the __spcosa__-package.
* The documentation of all functions has been fine tuned.
* A trick was included in the `optimMKV()`-function to avoid errors due to
  the LDLfactor error of the __gstat__-package.
* There also is a new function to compute the Pareto maximum and minimum values of the objective functions 
  that compose a multi-objective optimization problem (MOOP): `minmaxPareto()`.

## spsann 1.0-2.9000 (2015-10-27)

* Now `x.max` and `y.max` are, by default, set to half of the maximum distance
  in the x- and y-coordinates of `candi`, respectively. In the same manner, the
  argument `cutoff` of `optimPPL()` is set, by default, to half of the diagonal
  of the rectangle of sides `x.max` and `y.max`.

## spsann 1.0-2 (2015-10-25)

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

## spsann 1.0-1 (2015-07-30)

* Improved and updated documentation.
* __gstat__ is not a dependence any more.
* Fixed breaks due to changes in dependencies (__pedometrics__).

## spsann 1.0-0 (2015-07-14)

* Submission to CRAN.

## spsann 0.0.0.9012 (2015-07-14)

* An auxiliary function (`objSPSANN()`) was created to retrieve the energy state
  of an optimized sample configuration (OSC) at a given point of the 
  optimization.
* Long examples are not run any more to avoid overload of `R CMD check`.
* The authors' list was updated with the respective roles.

## spsann 0.0.0.9011 (2015-07-13)

* The documentation of all functions was significantly improved.
* Functions from default packages other than __base__ are now imported to
  comply with the new change to the CRAN policy described at
  http://developer.r-project.org/blosxom.cgi/R-devel/NEWS/2015/06/29#n2015-06-29.

## spsann 0.0.0.9010 (2015-07-05)

* Using `utils::globalVariables` to avoid the `R CMD check` note 
  `no visible binding for global variable [variable name]`. Source of the 
  solution: http://stackoverflow.com/a/12429344/3365410.
* The package __fields__ is not a dependency any more.
* New default values were attributed to the following arguments: `plotit`, 
  `track`, `verbose`, and `iteration`. The first three were set to `FALSE`, 
  while the last was set to `100`.
* `optimSPAN()` and `objSPAN()` are now full operational.

## spsann 0.0.0.9009 (2015-06-30)

* Several internal function were renamed using a pattern that includes the name
  of the respective objective function. For example, `.optimPPLcheck()` was 
  renamed as `.checkPPL()`, and `.getLagBreaks()` was renamed as `.lagsPPL()`.
  Note that the first part of the function name indicates what it does, while 
  the second indicates the objective function to which it applies. This 
  standardization is important to ease the construction of multi-objective
  optimization problems.

## spsann 0.0.0.9008 (2015-06-29)

* Improvements in the family of ACDC, CORR, and DIST functions.
* Several pairs of internal function that were originally designed to deal with
  different types of covariates (factor and numeric) were merged. Now a single
  function does the job by using the key argument `covars.type`.
* New internal functions now enable building multi-objective optimization
  problems more easily. They have also allowed to clean-up/simplify the source 
  code.
* A new `autofun` was created to set-up the covariates (`covar`).

## spsann 0.0.0.9007 (2015-06-12)

* The __rgeos__ and __plyr__ packages are not dependencies any more.
* The `boundary` of the spatial domain can now be estimated internally. The user 
  should use the __rgeos__ package if a more precise `boundary` is needed.
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

## spsann 0.0.0.9006 (2015-05-12)

* `spJitterFinite()` now tries to find an alternative point if the new point
  already is included in the sample. The number of tries is equal to the total
  number of points included in the sample. Because the more points we have, the
  more likely it is that the candidate point already is included in the sample.

## spsann 0.0.0.9005 (2015-04-29)

* `spJitterFinite()` now returns the old point if the new point already is in
  the sample. This is to avoid an infinite loop at the end of the optimization
  when the objective function creates a cluster of points.

## spsann 0.0.0.9004 (2015-04-24)

* New version of `optimACDC()`, including new argument definitions;
* In the multi-objective optimization problem case, now the graphical display 
  includes the many objective functions being optimized along with the utility
  function.

## spsann 0.0.0.9003 (2015-04-20)

* Special version designed for the course on Spatial Sampling for Mapping, 
  22 - 24 April 2015, Wageningen University Campus, Wageningen, The Netherlands,
  Under the auspices of the Graduate School for Production Ecology and Resource 
  Conservation (PE&RC).

## spsann 0.0.0.9002 (2015-04-19)

* new function to enable the user to define his/her own objective function;
* grammar check and enhanced documentation;

## spsann 0.0.0.9001 (2015-02-24)

* new functions derived from `optimACDC()`: `optimDIST()` and `optimCORR()`;
* new objective function: mean/maximum kriging variance;
* review of the family of PPL functions;
* using function tailored argument checking.

## spsann 0.0.0.9000 (2015-01-08)

* in-development package;
* importing functions from the package __pedometrics__;
* preparing documentation.
