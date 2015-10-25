## Changes
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

## Test environments
* local x86_64-pc-linux-gnu (ubuntu 14.04), R 3.2.2
* ubuntu 12.04 (on travis-ci), R 3.2.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies
There are no downstream dependencies of ***spsann***.
