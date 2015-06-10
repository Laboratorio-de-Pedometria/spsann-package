# Version 0.0.0.9007 (2015-06-10)
* BUG: `objACDC()` and `objCORR()` may not return the same criterion value
  of the optimized sample configuration if the number of iterations used
  in the optimization is equal to 100. The problem seems to disappear if a
  larger number of iterations is used.
* A directory called 'tools' was created, where R code chunks are included in
  individual files. These R code chunks are used in several functions of both
  families of `obj...()` and `optim...()` functions. They are used to 
  automatically build internal functions using `ASRtools::autoFunction()`.
  Currently, R code chunks are used to check the arguments of the family of 
  `optim...()` functions, prepare `points` and `candi`, set plotting options,
  estimate the `boundary`, and prepare for jittering.
* The `boundary` of the spatial domain can now be estimated intenally. As such,
  the `rgeos` package is not a dependency anymore. The user should relly on the
  `rgeos` package if a more precise `boundary` is required.
  

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
* importing functions from the package `pedometrics`;
* preparing documentation.
