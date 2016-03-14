## Changes
This is a major release of package ***spsann*** that includes several conceptual changes. Despite our efforts, 
it was not possible to guarantee the compatibility with previous versions. As such, the user will be warned 
about the changes upon attaching the package to the workspace. We have decided not to deprecate functions 
and function arguments because (1) this would require deprecating a lot of code, and (2) the user should 
first read the updated package documentation to understand the conceptual changes that we have made before
using it. This should not be a major problem because, as of now, the number of users of ***spsann*** appears to
be very limited -- all known users have been warned about this major release.

## Test environments
* local x86_64-pc-linux-gnu (ubuntu 14.04), R 3.2.3
* ubuntu 12.04.5 LTS (on travis-ci), R 3.2.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE.

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Alessandro Samuel-Rosa <alessandrosamuelrosa@gmail.com>'

Possibly mis-spelled words in DESCRIPTION:
  variogram (20:71)
```

## Downstream dependencies
There are no downstream dependencies of ***spsann***.
