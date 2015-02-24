#' spsann: optimization of sample patterns using spatial simulated annealing
#' 
#' The \strong{spsann} package is the tool for the optimization of sample
#' patterns using spatial simulated annealing in R. Spatial simulated annealing
#' is a well known method with widespread use to solve optimization problems in
#' the environmental sciences. This is mainly due to its robustness against 
#' local optima.
#' 
#' The \strong{spsann} functions offer many optimizing criteria for sampling for
#' variogram estimation (number of points or point-pairs per lag distance class 
#' - PPL), trend estimation (association/correlation and marginal distribution
#' of the covariates - ACDC), and spatial interpolation (mean squared shortest
#' distance - MSSD). \strong{spsann} also includes as an optimizing criterion 
#' the mean (or maximum) universal kriging variance (MUKV), which is useful when
#' the model of spatial variation is known. PPL, ACDC and MSSD were combined 
#' (PAN) for sampling when we are ignorant about the model of spatial variation
#' (\emph{terra incognita}).
#' 
#' \tabular{ll}{
#' Package: \tab spsann    \cr
#' Type:    \tab Package   \cr
#' Version: \tab 0.0.1.9000\cr
#' Date:    \tab 2015-02-24\cr
#' License: \tab GPL (>= 2)\cr
#' }
#' 
#' @name spsann-package
#' @aliases spsann-package spsann
#' @docType package
#' @author Author and Maintainer: Alessandro Samuel-Rosa
#' \email{alessandrosamuelrosa@@gmail.com}
#' 
#' @section Support:
#' The \strong{spsann} package was initially developed as part of the PhD
#' research project entitled \sQuote{Contribution to the Construction of Models
#' for Predicting Soil Properties}, developed by Alessandro Samuel-Rosa
#' under the supervision of Dr. Lúcia Helena Cunha dos Anjos (Universidade
#' Federal Rural do Rio de Janeiro, Brazil), Dr. Gustavo de Mattos Vasques
#' (Embrapa Solos, Brazil), and Dr. Gerard B. M. Heuvelink (ISRIC - World Soil
#' Information, the Netherlands). The project was/is supported from 2012 to 
#' 2016 by the CAPES Foundation, Ministry of Education of Brazil, and the CNPq
#' Foundation, Ministry of Science and Technology of Brazil.
#' 
#' @section Contributors:
#' Some of the solutions used to build the \strong{spsann} functions were found 
#' in the source code of other R-packages. The skeleton of the optimization
#' functions was adopted from the \pkg{intamapInteractive}, authored by
#' Edzer Pebesma <\email{edzer.pebesma@@uni-muenster.de}> and Jon Skoien 
#' <\email{jon.skoien@@gmail.com}>.
#' 
#' A few small solutions were adopted from the \pkg{SpatialTools} package,
#' authored by Joshua French <\email{joshua.french@@ucdenver.edu}>, and 
#' \pkg{clhs} package, authored by Pierre Roudier 
#' <\email{roudierp@@landcareresearch.co.nz}>.
#' 
#' Conceptual contributions were made by Gerard Heuvelink 
#' <\email{gerard.heuvelink@@wur.nl}> and Dick Brus <\email{dick.brus@@wur.nl}>.
#' 
#' @useDynLib spsann
NULL
