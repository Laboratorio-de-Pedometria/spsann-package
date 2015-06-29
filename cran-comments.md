## Test environments
* local x86_64-pc-linux-gnu, 3.2.1
* ubuntu 12.04 (on travis-ci), R 3.1.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking R code for possible problems ... NOTE
objACDC: no visible binding for global variable ‘covars.type’
objACDC: no visible binding for global variable ‘sm’
objACDC: no visible binding for global variable ‘n_pts’
objACDC: no visible binding for global variable ‘n_cov’
objACDC: no visible binding for global variable ‘n_candi’
objCORR: no visible binding for global variable ‘covars.type’
objCORR: no visible binding for global variable ‘sm’
objDIST: no visible binding for global variable ‘n_pts’
objDIST: no visible binding for global variable ‘covars.type’
objDIST: no visible binding for global variable ‘sm’
objDIST: no visible binding for global variable ‘n_cov’
objMKV: no visible binding for global variable ‘n_pts’
objPPL: no visible binding for global variable ‘n_pts’
optimACDC: no visible binding for global variable ‘covars.type’
optimACDC: no visible binding for global variable ‘sm’
optimACDC: no visible binding for global variable ‘n_pts’
optimACDC: no visible binding for global variable ‘n_cov’
optimACDC: no visible binding for global variable ‘n_candi’
optimACDC: no visible binding for global variable ‘wp’
optimCORR: no visible binding for global variable ‘covars.type’
optimCORR: no visible binding for global variable ‘sm’
optimCORR: no visible binding for global variable ‘wp’
optimDIST: no visible binding for global variable ‘n_pts’
optimDIST: no visible binding for global variable ‘covars.type’
optimDIST: no visible binding for global variable ‘sm’
optimDIST: no visible binding for global variable ‘n_cov’
optimDIST: no visible binding for global variable ‘wp’
optimMKV: no visible binding for global variable ‘n_pts’
optimMKV: no visible binding for global variable ‘wp’
optimMSSD: no visible binding for global variable ‘conf0’
optimMSSD: no visible binding for global variable ‘wp’
optimPPL: no visible binding for global variable ‘conf0’
optimPPL: no visible binding for global variable ‘n_pts’
optimPPL: no visible binding for global variable ‘wp’
optimSPAN: no visible binding for global variable ‘conf0’
optimSPAN: no visible binding for global variable ‘n_pts’
optimSPAN: no visible global function definition for ‘.covarsACDC’
optimSPAN: no visible global function definition for ‘.numStrata’
optimSPAN: no visible binding for global variable ‘n_candi’
optimSPAN: no visible binding for global variable ‘pre.distri’
optimSPAN: no visible global function definition for ‘.objNum’
optimSPAN : <anonymous>: no visible binding for global variable ‘n_candi’
optimSPAN: no visible binding for global variable ‘pop.prop’
optimSPAN: no visible global function definition for ‘.objFac’
optimSPAN: no visible binding for global variable ‘PAN’
optimSPAN: no visible binding for global variable ‘wp’
optimSPAN: no visible binding for global variable ‘ACDC’

These global variables are created by internal functions defined as expressions
which return more than on object.

## Downstream dependencies
There are no downstream dependencies of `spsann`
