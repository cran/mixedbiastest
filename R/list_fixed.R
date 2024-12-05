#' List Fixed Effects from an lmerMod Object
#'
#' This function lists the fixed effects coefficients from an `lmerMod` object, providing the index and name of each coefficient. This can help users construct contrast vectors (`k_list`) for use with the `mixedbiastest` function.
#'
#' @param model An object of class `lmerMod` fitted using `lmer` from the `lme4` package.
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{\code{Index}}{The index of each fixed effect coefficient.}
#'   \item{\code{Coefficient}}{The name of each fixed effect coefficient.}
#' }
#'
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview, which helped in constructing the structure of much of the code and the roxygen documentation. The code is based on the R code provided by Karl and Zimmerman (2020).
#'
#' @examples
#' if (requireNamespace("plm", quietly = TRUE) && requireNamespace("lme4", quietly = TRUE)) {
#'   library(lme4)
#'   data("Gasoline", package = "plm")
#'   # Fit a random effects model using lme4
#'   mixed_model <- lmer(lgaspcar ~ lincomep + lrpmg + lcarpcap + (1 | country),
#'                       data = Gasoline, REML = FALSE)
#'
#'   # List fixed effects
#'   fixed_effects <- list_fixed_effects(mixed_model)
#'   print(fixed_effects)
#'
#'   # Suppose we want to test the contrast: lincomep - lcarpvap
#'   p <- nrow(fixed_effects)
#'   k <- rep(0, p)
#'   k[fixed_effects$Index[fixed_effects$Coefficient == "lincomep"]] <- 1
#'   k[fixed_effects$Index[fixed_effects$Coefficient == "lcarpvap"]] <- -1
#'
#'   # Run the bias test with the custom contrast
#'   result <- mixedbiastest(mixed_model, k_list = list(k))
#'   print(result)
#'   plot(result)
#' } else {
#'   message("Please install 'plm' and 'lme4' packages to run this example.")
#' }

#' @importFrom lme4 fixef
#' @export
list_fixed_effects <- function(model) {
  if (!inherits(model, "lmerMod")) {
    stop("The 'model' must be an object of class 'lmerMod' from the 'lme4' package.")
  }
  coef_names <- names(fixef(model))
  data.frame(
    Index = seq_along(coef_names),
    Coefficient = coef_names,
    stringsAsFactors = FALSE
  )
}
