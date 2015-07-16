#  Template documentation for the family of CORR objective functions
################################################################################
#' @section Association/Correlation between covariates:
#' The \emph{correlation} between two numeric covariates is measured using the 
#' Pearson's r, a descriptive statistic that ranges from $-1$ to $+1$. 
#' This statistic is also known as the linear correlation coefficient.
#' 
#' When the set of covariates includes factor covariates, all numeric covariates 
#' are transformed into factor covariates. The factor levels are defined 
#' using the marginal sampling strata created from one of the two methods 
#' available (equal-area or equal-range strata).
#' 
#' The \emph{association} between two factor covariates is measured using the 
#' Cramér's v, a descriptive statistic that ranges from $0$ to $+1$. The closer 
#' to $+1$ the Cramér's v is, the stronger the association between two factor 
#' covariates. The main weakness of using the Cramér's v is that, while the 
#' Pearson's r shows the degree and direction of the association between two 
#' covariates (negative or positive), the Cramér's v only measures the degree 
#' of association (weak or strong).

