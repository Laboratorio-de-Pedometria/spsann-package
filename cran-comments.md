## Changes
Now ***spsann*** can be used to augment an existing sample configuration, that is, add new sampling points
to a spatial sample configuration generated using ***spsann*** or any other means. To do so, when using one 
of the functions from the family of `optim...()` functions, the user must pass to the function argument 
`points` an object of class `list` containing two named sub-arguments: `fixed`, a matrix with the 
coordinates of the existing sample configuration -- kept fixed during the optimization --, and `free`,
the number of sample points that should be added to the existing sample configuration -- free to move around
during the optimization.

## Test environments
* local x86_64-pc-linux-gnu (ubuntu 16.04), R 3.4.0
* ubuntu 12.04.5 LTS (on travis-ci), R 3.4.0
* win-builder -- R Under development (unstable) (2017-06-23 r72844)

## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE.

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Alessandro Samuel-Rosa <alessandrosamuelrosa@gmail.com>'

Possibly mis-spelled words in DESCRIPTION:
  multi (23:5)
  variogram (20:61)
```

## Downstream dependencies
There are no downstream dependencies of ***spsann***.
