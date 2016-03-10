#' The spsann package
#' 
#' Optimization of sample configurations using spatial simulated annealing
#' 
#' @details 
#' \pkg{spsann} is a package for the optimization of spatial sample configurations using spatial simulated 
#' annealing. It includes multiple functions with different objective functions to optimize spatial sample 
#' configurations for various purposes such as variogram estimation, spatial trend estimation, and spatial 
#' interpolation. Most of the objective functions were designed to optimize spatial sample configurations when
#' we know very little about the model of spatial variation and sampling is limited to a single phase.
#' 
#' A general purpose function enables to user to define his/her own objective function and plug it in the
#' spatial simulated annealing algorithm.
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
#' Version: \tab 1.0-2.9010\cr
#' Date:    \tab 2015-03-10\cr
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
#' \pkg{spsann} was initially developed as part of the PhD research project entitled \sQuote{Contribution to 
#' the Construction of Models for Predicting Soil Properties}, developed by Alessandro Samuel-Rosa under the 
#' supervision of Lúcia Helena Cunha dos Anjos (Universidade Federal Rural do Rio de Janeiro, Brazil), Gustavo 
#' de Mattos Vasques (Embrapa Solos, Brazil), and Gerard B. M. Heuvelink (ISRIC -- World Soil Information, the 
#' Netherlands). The project was supported from March/2012 to February/2016 by the CAPES Foundation, Ministry 
#' of Education of Brazil, and the CNPq Foundation, Ministry of Science and Technology of Brazil.
#' 
#' @section Contributors:
#' Some of the solutions used to build \pkg{spsann} were found in the source code of other R-packages and 
#' scripts developed and published by other researchers. For example, the original skeleton of the 
#' optimization functions was adopted from the \pkg{intamapInteractive} package with the approval of the 
#' package authors, Edzer Pebesma \email{edzer.pebesma@@uni-muenster.de} and Jon Skoien 
#' \email{jon.skoien@@gmail.com}. The current skeleton is based on the later adoption of several solutions 
#' implemented in the script developed and published by Murray Lark \email{mlark@@bgs.ac.uk} as part of a
#' short course offered for the first time at the 2015 EGU General Assembly in Vienna, Austria.
#' 
#' A few small solutions were adopted from the packages \pkg{SpatialTools}, authored by Joshua French 
#' \email{joshua.french@@ucdenver.edu}, \pkg{clhs}, authored by Pierre Roudier 
#' \email{roudierp@@landcareresearch.co.nz}, and \pkg{spcosa}, authored by Dennis Walvoort 
#' \email{dennis.Walvoort@@wur.nl}, Dick Brus \email{dick.brus@@wur.nl}, and Jaap de Gruijter 
#' \email{Jaap.degruijter@@wur.nl}.
#' 
#' Major conceptual contributions were made by Gerard Heuvelink \email{gerard.heuvelink@@wur.nl}, Dick Brus 
#' \email{dick.brus@@wur.nl}, Murray Lark \email{mlark@@bgs.ac.uk}, and Edzer Pebesma 
#' \email{edzer.pebesma@@uni-muenster.de}.
#' 
#' @useDynLib spsann
NULL
