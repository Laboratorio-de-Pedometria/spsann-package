#' spsann: Optimization of Sample Configurations using Spatial Simulated
#' Annealing
#' 
#' Methods to optimize sample configurations using spatial simulated annealing. 
#' Multiple objective functions are implemented for various purposes, such as 
#' variogram estimation, trend estimation, and spatial interpolation. A general 
#' purpose spatial simulated annealing function enables the user to define 
#' his/her own objective function.
#' 
#' \pkg{spsann} is the R package for the optimization of sample configurations 
#' using spatial simulated annealing. It includes multiple functions with 
#' different objective functions to optimize sample configurations for 
#' variogram estimation (number of points or point-pairs per lag distance 
#' class), trend estimation (association/correlation and marginal distribution
#' of the covariates), and spatial interpolation (mean squared shortest 
#' distance). \strong{spsann} also includes objective functions that can be 
#' used when the model of spatial variation is known (mean (or maximum) kriging
#' variance). Some objective functions were combined to optimize sample 
#' configurations when we are ignorant about the model of spatial variation 
#' (\emph{terra incognita}). A general purpose function enables to user to 
#' define his/her own objective function and plug it in the spatial simulated 
#' annealing algorithm.
#' 
#' Spatial simulated annealing is a well known method with widespread use to 
#' solve optimization problems in the environmental sciences. This is mainly 
#' due to its robustness against local optima. At each iteration, the algorithm
#' evaluates if a worsening solution can be accepted. The chance of accepting
#' worsening solutions reduces as the number of iterations increases. 
#' 
#' \tabular{ll}{
#' Package: \tab spsann    \cr
#' Type:    \tab Package   \cr
#' Version: \tab 0.0.0.9007\cr
#' Date:    \tab 2015-06-10\cr
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
#' \pkg{spsann} was initially developed as part of the PhD
#' research project entitled \sQuote{Contribution to the Construction of Models
#' for Predicting Soil Properties}, developed by Alessandro Samuel-Rosa
#' under the supervision of Lúcia Helena Cunha dos Anjos (Universidade
#' Federal Rural do Rio de Janeiro, Brazil), Gustavo de Mattos Vasques
#' (Embrapa Solos, Brazil), and Gerard B. M. Heuvelink (ISRIC - World Soil
#' Information, the Netherlands). The project was/is supported from 2012 to 
#' 2016 by the CAPES Foundation, Ministry of Education of Brazil, and the CNPq
#' Foundation, Ministry of Science and Technology of Brazil.
#' 
#' @section Contributors:
#' Some of the solutions used to build \strong{spsann} were found 
#' in the source code of other R-packages. The skeleton of the optimization
#' functions was adopted from the \pkg{intamapInteractive} package, authored by
#' Edzer Pebesma \email{edzer.pebesma@@uni-muenster.de} and Jon Skoien 
#' \email{jon.skoien@@gmail.com}.
#' 
#' A few small solutions were adopted from the \pkg{SpatialTools} package,
#' authored by Joshua French \email{joshua.french@@ucdenver.edu}, and 
#' \pkg{clhs} package, authored by Pierre Roudier 
#' \email{roudierp@@landcareresearch.co.nz}.
#' 
#' Conceptual contributions were made by Gerard Heuvelink 
#' \email{gerard.heuvelink@@wur.nl}, Dick Brus \email{dick.brus@@wur.nl},
#' Murray Lark \email{mlark@@bgs.ac.uk}, and Edzer Pebesma 
#' \email{edzer.pebesma@@uni-muenster.de}.
#' 
#' @useDynLib spsann
NULL
