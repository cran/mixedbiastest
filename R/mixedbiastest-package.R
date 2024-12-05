#' mixedbiastest: Bias Diagnostics for Linear Mixed Models
#'
#' The `mixedbiastest` package provides a function to perform bias diagnostics on linear mixed models fitted with `lmer` from the `lme4` package. It implements permutation tests for assessing the bias of fixed effects, as described in Karl and Zimmerman (2021).
#'
#' @section Functions:
#' \describe{
#'   \item{\code{\link{mixedbiastest}}}{Performs the bias diagnostic test.}
#'   \item{\code{\link{print.mixedbiastest}}}{Prints the results of the bias diagnostic.}
#'   \item{\code{\link{plot.mixedbiastest}}}{Plots the permutation distributions and observed test statistics for each fixed effect.}
#'   \item{\code{\link{list_fixed_effects}}}{ List Fixed Effects from an lmerMod Object.}
#' }
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview, which helped in constructing the structure of much of the code and the roxygen documentation. The code is based on the R code provided by Karl and Zimmerman (2020).
#'
#' @references
#' Karl, A. T., & Zimmerman, D. L. (2021). A diagnostic for bias in linear mixed model estimators induced by dependence between the random effects and the corresponding model matrix. \emph{Journal of Statistical Planning and Inference}, \emph{212}, 70–80. \doi{10.1016/j.jspi.2020.06.004}
#'
#' Karl, A., & Zimmerman, D. (2020). Data and Code Supplement for 'A Diagnostic for Bias in Linear Mixed Model Estimators Induced by Dependence Between the Random Effects and the Corresponding Model Matrix'. Mendeley Data, V1. \doi{10.17632/tmynggddfm.1}
#'
#' @docType package
#' @name mixedbiastest-package
"_PACKAGE"
