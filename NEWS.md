# Version 1.0.1 (2015-07-30)
* Improved and updated documentation.
* ***gstat*** is not a dependence any more.
* Fixed breaks due to changes in dependencies (***pedometrics***).

# Version 1.0.0 (2015-07-14)
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
