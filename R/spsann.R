#' Optimization of sample patterns using spatial simulated annealing
#' 
#' \bold{spsann} offers many optimizing criteria for sampling for variogram
#' estimation (PPL), trend estimation (ACDC), and spatial interpolation (MSSD).
#' \bold{spsann} also includes an optimizing criterion for when the model of
#' spatial variation is known (MUKV). PPL, ACDC and MSSD were combined (PAN) for
#' sampling when we are ignorant about the model of spatial variation.
#' 
#' \tabular{ll}{
#' Package: \tab spsann    \cr
#' Type:    \tab Package   \cr
#' Version: \tab 1.0       \cr
#' Date:    \tab 2015-01-07\cr
#' License: \tab GPL (>= 2)\cr
#' }
#' 
#' @name spsann-package
#' @aliases spsann-package spsann
#' @docType package
#' @author Author and Maintainer: Alessandro Samuel-Rosa
#' \email{alessandrosamuelrosa@@gmail.com}
#' 
#' Thesis advisors: Lúcia Helena Cunha dos Anjos, Gustavo de Mattos Vasques,
#' Gerard B. M. Heuvelink
#' 
#' Contributors: Edzer Pebesma, Jon Skoien, Joshua French, Ken Kleinman, 
#' Dick Brus
#' 
#' \bold{spsann} was developed as part of the PhD research project
#' (2012-2016) entitled Contribution to the Construction of Models for 
#' Predicting Soil Properties, developed by Alessandro Samuel-Rosa under the
#' supervision of Dr Lúcia Helena Cunha dos Anjos (Universidade Federal Rural do
#' Rio de Janeiro, Brazil), Dr Gustavo de Mattos Vasques (Embrapa Solos, 
#' Brazil), and Dr Gerard B. M. Heuvelink (ISRIC - World Soil Information, the
#' Netherlands). The project is supported by CAPES (Minitério da Educação) and
#' CNPq (Ministério da Ciência, Tecnologia e Inovação).
#' @useDynLib spsann
NULL
