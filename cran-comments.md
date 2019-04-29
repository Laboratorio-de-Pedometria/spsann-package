## Changes

The new version of the __spsann__ package includes some bug fixes and a few modifications. Users now can 
choose how `optimCLHS` computes objective function values: as in the original paper or as in the FORTRAN 
implementation. Users now also must inform the `weights` passed to `optimCLHS` as to guarantee that s/he is 
aware of what s/he is doing. The same apples to other functions that deal with multi-objective optimization 
problems: `optimACDC` and `optimSPAM`. Another important modification in the current version of __spsann__ is 
the possibility to use a finite set of candidate locations by setting `cellsize = 0`. This is useful when 
optimizing sample points only in the feature space and should reduce the computation time needed to find the
solution.

## Test environments

* OK: local, x86_64-pc-linux-gnu (64-bit), Ubuntu 18.04.01 LTS, R 3.6.0
* OK: travis-ci, x86_64-pc-linux-gnu (64-bit), Ubuntu 14.04.5 LTS, R 3.5.3
* OK: rhub, Debian Linux, R-devel, GCC ASAN/UBSAN
* OK: rhub, Fedora Linux, R-devel, clang, gfortran
* WARNING: rhub, Ubuntu Linux 16.04 LTS, R-release, GCC
* NOTE: rhub, Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* OK: winbuilder, x86_64-w64-mingw32 (64-bit), Windows, R 3.5.3
* OK: winbuilder, x86_64-w64-mingw32 (64-bit), Windows, R 3.6.0
* OK: winbuilder, x86_64-w64-mingw32 (64-bit), Windows, R 3.6.0 RC

## R CMD check results

There were no ERRORs.

There were two NOTEs and one WARNING.

### Notes

On rhub, Windows Server 2008 R2 SP1, R-devel, 32/64 bit:

```
* checking sizes of PDF files under 'inst/doc' ... NOTE
Unable to find GhostScript executable to run checks on size reduction
```

This note appears to be due to the `rhub` configuration and issues have been registered before for different
platforms:

* https://github.com/r-hub/rhub/issues/82
* https://github.com/r-hub/rhub/issues/51

A new issue was registered at https://github.com/r-hub/rhub/issues/263 and the `rhub` package maintainer will
look into it.

On all platforms:

```
─  checking CRAN incoming feasibility ... Note_to_CRAN_maintainers (2.2s)
   Maintainer: 'Alessandro Samuel-Rosa <alessandrosamuelrosa@gmail.com>'
```

This NOTE can be ignored.

### Warnings

On rhub, Ubuntu Linux 16.04 LTS, R-release, GCC:

```
* checking compilation flags used ... WARNING
Compilation used the following non-portable flag(s):
  ‘-Wdate-time’ ‘-Werror=format-security’ ‘-Wformat’
```

This WARNING can be ignored as stated at https://stat.ethz.ch/pipermail/r-package-devel/2018q3/002889.html.

## Downstream dependencies

There are no downstream dependencies of ***spsann***.
