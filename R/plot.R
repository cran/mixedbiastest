#' Plot Method for Bias Diagnostic Results
#'
#' Plots the permutation distributions and observed test statistics for each fixed effect.
#'
#' @param x An object of class `"mixedbiastest"`.
#' @param ... Additional arguments (currently not used).
#' @return A \code{ggplot} object showing permutation distributions for all fixed effects.
#' @method plot mixedbiastest
#' @importFrom ggplot2 ggplot geom_histogram geom_vline facet_wrap labs theme_minimal geom_text aes
#' @importFrom stats rnorm
#' @importFrom rlang .data
#' @export
plot.mixedbiastest <- function(x, ...) {
  permutation_values_list <- x$permutation_values
  results <- x$results_table
  k_list_names <- results$Fixed_Effect

  # Combine all permutation values into a single data frame
  df_list <- list()
  annotations <- data.frame(Fixed_Effect = character(), label = character(), stringsAsFactors = FALSE)

  for (fe_name in k_list_names) {
    permutation_values <- permutation_values_list[[fe_name]]
    observed_value <- as.numeric(results$Bias_Estimate[results$Fixed_Effect == fe_name])

    # Check if all permutation values are the same
    if (length(unique(permutation_values)) == 1) {
      # Indicate that permutation values are constant
      annotations <- rbind(annotations, data.frame(
        Fixed_Effect = fe_name,
        label = "No variation in permutation values"
      ))
      permutation_values_current <- permutation_values
    } else {
      permutation_values_current <- permutation_values
    }

    df <- data.frame(
      PermutationValues = permutation_values_current,
      Fixed_Effect = fe_name,
      Observed_Value = observed_value
    )
    df_list[[fe_name]] <- df
  }

  combined_df <- do.call(rbind, df_list)

  # Adjust the factor levels
  combined_df$Fixed_Effect <- factor(combined_df$Fixed_Effect, levels = k_list_names)
  annotations$Fixed_Effect <- factor(annotations$Fixed_Effect, levels = k_list_names)



  # Plot using facets
  p <- ggplot(combined_df, aes(x = .data$PermutationValues)) +
    geom_histogram( fill = "lightblue", color = "black", na.rm = TRUE) +
    geom_vline(aes(xintercept = .data$Observed_Value), color = "red",
               linetype = "dashed", linewidth = 1) +
    facet_wrap(~ Fixed_Effect, scales = "free") +
    labs(title = "Bias Diagnostic for Fixed Effects",
         x = expression(hat(nu) %*% hat(eta)),
         y = "Frequency") +
    theme_minimal()

  # Add annotations if any
  if (nrow(annotations) > 0) {
    p <- p + geom_text(
      data = annotations,
      aes(
        x = mean(combined_df$PermutationValues, na.rm = TRUE),
        y = 0,
        label = .data$label
      ),
      color = "red",
      size = 3,
      hjust = 0.5,
      vjust = -1,
      inherit.aes = FALSE
    )
  }

  # Suppress warnings and print the plot
  suppressWarnings(print(p))

  # Optionally return the plot if desired
  # return(p)
}
