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

* new functions derived from `optimACDC`: `optimDIST` and `optimCORR`;
* new objective function: mean/maximum kriging variance;
* review of the family of PPL functions;
* using function tailored argument checking.

# Version 0.0.0.9000 (2015-01-08)

* in-development package;
* importing functions from the package `pedometrics`;
* preparing documentation.
