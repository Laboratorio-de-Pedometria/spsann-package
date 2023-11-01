# Changes

This is a patch. Only a few tests were performed on various platforms. A single change was made in
the vignette file, because the option 'english' passed to the babel package was resulting in an
error in the vignette build process.

## Test environments

* OK: local, x86_64-pc-linux-gnu (64-bit), Ubuntu 22.04.3 LTS, R 4.3.1
* OK: winbuilder, x86_64-w64-mingw32 (64-bit), Windows, R 4.3.2
* OK: winbuilder, x86_64-w64-mingw32 (64-bit), Windows, R 4.2.3
* OK: winbuilder, x86_64-w64-mingw32 (64-bit), Windows, R Under development (unstable)
* WARNING: rhub, Fedora Linux, R-devel, clang, gfortran
* ERROR: rhub, Windows Server 2022, R-devel, 64 bit
* PREPERROR: rhub, Debian Linux, R-devel, clang, ISO-8859-15 locale
* PREPERROR: rhub, Ubuntu Linux 20.04.1 LTS, R-release, GCC

## R CMD check results

### Errors

There was an ERROR in the `rhub` platform, Windows Server 2022, R-devel, 64 bit:

```
* checking re-building of vignette outputs ... [28s] ERROR
Error(s) in re-building vignettes:
  ...
--- re-building 'spsann.Rnw' using knitr
Warning in system(paste(shQuote(texi2dvi), if (quiet) "--quiet" else "",  :
  running command '"C:\PROGRA~1\MiKTeX\miktex\bin\x64\texify.exe" --quiet --pdf "spsann.tex" --max-iterations=20 -I "C:/PROGRA~1/R/R-devel/share/texmf/tex/latex" -I "C:/PROGRA~1/R/R-devel/share/texmf/bibtex/bst"' had status 1
Error: processing vignette 'spsann.Rnw' failed with diagnostics:
Failed to locate 'texi2pdf' output file 'spsann.pdf' for vignette with name 'spsann' and engine 'knitr::knitr'. The following files exist in working directory 'C:\Users\USERDjNWOvwxAt\spsann.Rcheck\vign_test\spsann\vignettes': 'spsann-concordance.tex' (133 bytes), 'spsann.R' (6150 bytes), 'spsann.Rnw' (15438 bytes), 'spsann.fdb_latexmk' (13135 bytes), 'spsann.fls' (42324 bytes), 'spsann.log' (4086 bytes), 'spsann.tex' (24036 bytes)
--- failed re-building 'spsann.Rnw'

SUMMARY: processing the following file failed:
  'spsann.Rnw'

Error: Vignette re-building failed.
Execution halted
```

This error suggests that it is necessary to check the vignette and how it is built.

### Warnings

There was one warning on rhub, Fedora Linux, R-devel, clang, gfortran:

```
* checking re-building of vignette outputs ... WARNING
Error(s) in re-building vignettes:
  ...
--- re-building ‘spsann.Rnw’ using knitr
Error: processing vignette 'spsann.Rnw' failed with diagnostics:
Running 'texi2dvi' on 'spsann.tex' failed.
LaTeX errors:
! LaTeX Error: File `framed.sty' not found.

Type X to quit or <RETURN> to proceed,
or enter new name. (Default extension: sty)

! Emergency stop.
<read *> 
         
l.27 \makeatletter
                  ^^M
!  ==> Fatal error occurred, no output PDF file produced!
--- failed re-building ‘spsann.Rnw’

SUMMARY: processing the following file failed:
  ‘spsann.Rnw’

Error: Vignette re-building failed.
Execution halted
```

This warning suggests that it is necessary to check the vignette and how it is built.

### Notes

The following notes appeared on virtually all platforms:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Alessandro Samuel-Rosa <alessandrosamuelrosa@gmail.com>'

New submission

Package was archived on CRAN

Possibly misspelled words in DESCRIPTION:
  variogram (35:61)

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2021-10-06 as requires archived package
    'pedometrics'.

Found the following (possibly) invalid URLs:
  URL: https://cran.r-project.org/web/checks/check_results_spsann.html
    From: README.md
    Status: 404
    Message: Not Found

Found the following (possibly) invalid file URI:
  URI: commits/master
    From: README.md

checking Rd cross-references ... NOTE
  Packages unavailable to check Rd xrefs: ‘raster’, ‘geoR’, ‘spcosa’
```

These NOTEs can be ignored for now.

## Downstream dependencies

There are no downstream dependencies of ***spsann***.
