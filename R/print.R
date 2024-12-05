#' Print Method for Bias Diagnostic Results
#'
#' Prints the results of the bias diagnostic in a formatted table.
#'
#' @param x An object of class `"mixedbiastest"`.
#' @param ... Additional arguments (currently not used).
#' @return No return value. This function is called for its side effects (printing the results).
#' @method print mixedbiastest
#' @export
print.mixedbiastest <- function(x, ...) {
  cat("Bias Diagnostic Results:\n")
  print(x$results_table, row.names = FALSE)
}
