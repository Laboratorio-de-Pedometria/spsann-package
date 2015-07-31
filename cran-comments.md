## Changes
* Fixed breaks due to changes in dependencies (***pedometrics***).

## Test environments
* local x86_64-pc-linux-gnu (ubuntu 14.04), R 3.2.1
* ubuntu 12.04 (on travis-ci), R 3.2.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There were two NOTEs:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Alessandro Samuel-Rosa <alessandrosamuelrosa@gmail.com>'

Possibly mis-spelled words in DESCRIPTION:
  variogram (31:23)
```

```
* checking package dependencies ... NOTE
  No repository set, so cyclic dependency check skipped
```

These notes do not affect the behaviour of the package in any of the tested
environments.

## Downstream dependencies
There are no downstream dependencies of ***spsann***.
