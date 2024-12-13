#' Bias Diagnostic for Linear Mixed Models
#'
#' Performs a permutation test to assess the bias of fixed effects in a linear mixed model fitted with `lmer`.
#' This function computes the test statistic and performs the permutation test, returning an object of class `"mixedbiastest"`.
#'
#' @param model An object of class `lmerMod` fitted using `lmer` from the `lme4` package.
#' @param n_permutations Integer. Number of permutations to perform (default is 10000).
#' @param k_list Optional list of numeric vectors. Each vector specifies a linear combination of fixed effects to test. If `NULL`, each fixed effect is tested individually.
#' @param verbose Logical. If `TRUE`, prints detailed messages during execution.
#' @details
#' **Note:** This function currently supports only models with diagonal random effects covariance matrices (i.e., the G matrix is diagonal). The methodology for non-diagonal G matrices is described in Karl and Zimmerman (2021), but is not implemented in this version of the package.
#'
#' See the \code{\link{list_fixed_effects}} function if you would like to estimate the bias of a contrast of fixed effects.
#'
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview, which helped in constructing the structure of much of the code and the roxygen documentation. The code is based on the R code provided by Karl and Zimmerman (2020).
#'
#' @return An object of class \code{"mixedbiastest"} containing:
#'
#' \describe{
#'   \item{\code{results_table}}{A data frame with the test results for each fixed effect or contrast, including bias estimates and p-values.}
#'   \item{\code{permutation_values}}{A list of numeric vectors containing permutation values for each fixed effect or contrast.}
#'   \item{\code{model}}{The original \code{lmerMod} model object provided as input.}
#' }
#' @references
#' Karl, A. T., & Zimmerman, D. L. (2021). A diagnostic for bias in linear mixed model estimators induced by dependence between the random effects and the corresponding model matrix. \emph{Journal of Statistical Planning and Inference}, \emph{212}, 70–80. \doi{10.1016/j.jspi.2020.06.004}
#'
#' Karl, A., & Zimmerman, D. (2020). Data and Code Supplement for 'A Diagnostic for Bias in Linear Mixed Model Estimators Induced by Dependence Between the Random Effects and the Corresponding Model Matrix'. Mendeley Data, V1. \doi{10.17632/tmynggddfm.1}
#'
#' @examples
#' if (requireNamespace("plm", quietly = TRUE) && requireNamespace("lme4", quietly = TRUE)) {
#' library(lme4)
#' data("Gasoline", package = "plm")
#' # Fit a random effects model using lme4
#' mixed_model <- lmer(lgaspcar ~ lincomep + lrpmg + lcarpcap + (1 | country),
#'                     data = Gasoline)
#' result <- mixedbiastest(mixed_model)
#' print(result)
#' plot(result)
#' }
#' if (requireNamespace("lme4", quietly = TRUE)) {
#' library(lme4)
#' example_model <- lmer(y ~ x + (1| group), data = example_data)
#' result2 <- mixedbiastest(example_model)
#' print(result2)
#' plot(result2)
#'
#'
#' #Simulate data
#' set.seed(123)
#' n_groups <- 30
#' n_obs_per_group <- 10
#' group <- rep(1:n_groups, each = n_obs_per_group)
#' x <- runif(n_groups * n_obs_per_group)
#' beta0 <- 2
#' beta1 <- 5
#' sigma_u <- 1
#' sigma_e <- 0.5
#' u <- rnorm(n_groups, 0, sigma_u)
#' e <- rnorm(n_groups * n_obs_per_group, 0, sigma_e)
#' y <- beta0 + beta1 * x + u[group] + e

#' data_sim <- data.frame(y = y, x = x, group = factor(group))

# Fit an lmer model (random intercept)
#' model3 <- lmer(y ~ x + (1 | group), data = data_sim)

# Perform the bias diagnostic with verbose output
#' result3 <- mixedbiastest(model3, verbose = TRUE)
#' plot(result3)
#'}
#' @importFrom lme4 getME fixef VarCorr
#' @importFrom stats sigma
#' @import Matrix
#' @export
mixedbiastest <- function(model, n_permutations = 10000, k_list = NULL, verbose = FALSE) {

  # Check if model is of class lmerMod
  if (!inherits(model, "lmerMod")) {
    stop("The 'model' must be an object of class 'lmerMod' from the 'lme4' package.")
  }

  # Check for singular fit
  if (lme4::isSingular(model, tol = 1e-4)) {
    warning("The model fit is singular. Some variance components are zero.")
  }

  # Extract model components
  X <- lme4::getME(model, "X")       # Fixed effects design matrix
  Z <- lme4::getME(model, "Z")       # Random effects design matrix
  y <- lme4::getME(model, "y")       # Response vector
  beta_hat <- lme4::fixef(model)     # Estimated fixed effects
  b_hat <- lme4::getME(model, "b")   # Estimated random effects (as a vector)

  # Residual variance estimate
  sigma2_hat <- sigma(model)^2  # Use residual variance estimate

  # Extract random effects terms and grouping factors
  VarCorr_list <- lme4::VarCorr(model)  # Variance-covariance matrices of random effects
  cnms <- lme4::getME(model, "cnms")     # Names of random effects
  flist <- lme4::getME(model, "flist")   # Grouping factors
  Gp <- lme4::getME(model, "Gp")         # Group pointers into b_hat

  if (verbose) {
    message("Dimensions of X (fixed effects design matrix): ", paste(dim(X), collapse = " x "))
    message("Dimensions of Z (random effects design matrix): ", paste(dim(Z), collapse = " x "))
    message("Class of Z: ", class(Z))
    message("Length of y (response vector): ", length(y))
    message("Length of beta_hat (estimated fixed effects): ", length(beta_hat))
    message("Length of b_hat (estimated random effects): ", length(b_hat))
    message("Estimated residual variance (sigma^2): ", sigma2_hat)
  }

  # Number of observations and parameters
  n <- length(y)
  p <- length(beta_hat)

  # Initialize G and G_inv as block-diagonal sparse matrices
  G_blocks <- list()
  G_inv_blocks <- list()
  b_hat_list <- list()

  # Helper function to check if a matrix is diagonal
  isDiagonal <- function(mat) {
    return(all(abs(mat[lower.tri(mat)]) < .Machine$double.eps & abs(mat[upper.tri(mat)]) < .Machine$double.eps))
  }

  # Loop over random effects terms
  for (i in seq_along(cnms)) {
    # The grouping factor for the i-th random effects term
    grouping_factor_name <- names(cnms)[i]
    n_re_per_group <- length(cnms[[i]])
    n_groups <- nlevels(flist[[grouping_factor_name]])

    # Indices for b_hat
    start_idx <- Gp[i] + 1
    end_idx <- Gp[i + 1]

    # Extract the portion of b_hat for this term
    b_hat_term <- b_hat[start_idx:end_idx]
    b_hat_matrix <- matrix(b_hat_term, nrow = n_groups, ncol = n_re_per_group, byrow = TRUE)
    b_hat_list[[i]] <- b_hat_matrix  # Use index as key

    # Extract the corresponding covariance matrix
    cov_mat <- as.matrix(VarCorr_list[[grouping_factor_name]])

    if (verbose) {
      message("Term: ", i)
      message("  Grouping factor: ", grouping_factor_name)
      message("  Covariance matrix dimensions: ", paste(dim(cov_mat), collapse = " x "))
      message("  Number of random effects per group: ", n_re_per_group)
      message("  Number of groups: ", n_groups)
    }

    # Skip terms with zero variance components
    if (all(diag(cov_mat) < .Machine$double.eps)) {
      if (verbose) {
        message("Skipping random effects term ", i, " due to zero variance.")
      }
      next
    }

    # Check if cov_mat is diagonal
    if (!isDiagonal(cov_mat)) {
      stop(paste("The covariance matrix for random effect term", i, "is not diagonal.",
                 "This function currently supports only models with diagonal random effects covariance matrices."))
    }

    # Build G blocks
    diag_elements <- diag(cov_mat)
    n_b <- n_re_per_group * n_groups
    G_block <- Diagonal(n = n_b, x = rep(diag_elements, times = n_groups))
    G_blocks[[i]] <- G_block

    # Build G_inv blocks
    G_inv_diag_elements <- 1 / diag_elements
    G_inv_block <- Diagonal(n = n_b, x = rep(G_inv_diag_elements, times = n_groups))
    G_inv_blocks[[i]] <- G_inv_block

    if (verbose) {
      message("  Dimensions of G_block: ", paste(dim(G_block), collapse = " x "))
      message("  Extracted b_hat indices: ", start_idx, " to ", end_idx)
      message("  Dimensions of b_hat_matrix: ", paste(dim(b_hat_matrix), collapse = " x "))
    }
  }

  # Combine G and G_inv blocks
  G <- bdiag(G_blocks)
  G_inv <- bdiag(G_inv_blocks)

  if (verbose) {
    message("Dimensions of combined G matrix: ", paste(dim(G), collapse = " x "))
  }

  # Construct R_inv
  R_inv <- (1 / sigma2_hat) * Diagonal(n)

  # Compute T = (Z' R_inv Z + G_inv)^{-1}
  Z_R_inv <- (1 / sigma2_hat) * Z
  T_mat <- crossprod(Z_R_inv, Z) + G_inv
  T_mat <- symmpart(T_mat)  # Ensure symmetry
  chol_T <- chol(T_mat)
  T <- chol2inv(chol_T)

  # Compute V_inv = R_inv - R_inv Z T Z' R_inv
  V_inv <- R_inv - Z_R_inv %*% T %*% t(Z_R_inv)

  if (verbose) {
    message("Dimensions of V_inv matrix: ", paste(dim(V_inv), collapse = " x "))
  }

  # Compute XVX_inv
  XVX <- crossprod(X, V_inv %*% X)
  XVX <- symmpart(XVX)  # Ensure symmetry (important before inversion)
  chol_XVX <- chol(XVX)
  XVX_inv <- chol2inv(chol_XVX)

  if (verbose) {
    message("Dimensions of XVX_inv matrix: ", paste(dim(XVX_inv), collapse = " x "))
  }

  # If k_list is NULL, create default k vectors for each fixed effect
  if (is.null(k_list)) {
    k_list <- lapply(1:p, function(j) {
      k <- rep(0, p)
      k[j] <- 1
      return(k)
    })
    names(k_list) <- names(beta_hat)
  } else {
    # Validate k_list
    if (!is.list(k_list)) {
      stop("k_list must be a list of numeric vectors.")
    }
    for (i in seq_along(k_list)) {
      k <- k_list[[i]]
      if (!is.numeric(k) || length(k) != p) {
        stop("Each k vector in k_list must be numeric and of length equal to the number of fixed effects.")
      }
    }
  }

  # Assign default names if names(k_list) is NULL
  if (is.null(names(k_list))) {
    names(k_list) <- paste0("Contrast_", seq_along(k_list))
  }

  # Initialize storage for results
  results <- data.frame(
    Fixed_Effect = names(k_list),
    #Observed_Value = NA,
    Mean_Permuted_Value = NA,
    P_Value = NA,
    Bias_Estimate = NA,
    stringsAsFactors = FALSE
  )

  # Initialize list to store permutation values
  permutation_values_list <- list()

  # Loop over each k vector
  for (j in seq_along(k_list)) {
    k <- k_list[[j]]
    fe_name <- names(k_list)[j]

    # Compute nu_hat = k' * (X' V^{-1} X)^{-1} * X' V^{-1} Z
    nu_hat <- as.numeric(k %*% XVX_inv %*% crossprod(X, V_inv %*% Z))

    # Compute the observed value: nu_hat' * b_hat
    observed_value <- sum(nu_hat * b_hat)

    if (verbose) {
      message("\nComputing bias diagnostic for fixed effect: ", fe_name)
      message("  Observed value (nu_hat' * b_hat): ", observed_value)
    }

    # Perform permutation test
    permutation_values <- numeric(n_permutations)

    for (i_perm in 1:n_permutations) {
      permuted_b_hat <- numeric(length(b_hat))
      for (k_term in seq_along(cnms)) {
        # Skip terms that were excluded due to zero variance
        if (length(b_hat_list) < k_term || is.null(b_hat_list[[k_term]])) {
          next
        }
        b_matrix <- b_hat_list[[k_term]]
        n_rows <- nrow(b_matrix)
        n_cols <- ncol(b_matrix)
        permuted_b_matrix <- b_matrix
        for (col_idx in 1:n_cols) {
          permuted_b_matrix[, col_idx] <- sample(b_matrix[, col_idx])
        }
        # Indices for permuted_b_hat
        idx_start <- Gp[k_term] + 1
        idx_end <- Gp[k_term + 1]
        permuted_b_hat[idx_start:idx_end] <- as.vector(t(permuted_b_matrix))
      }
      # Compute the test statistic
      permutation_values[i_perm] <- sum(nu_hat * permuted_b_hat)
    }

    # Ensure permutation_values are finite
    permutation_values <- permutation_values[is.finite(permutation_values)]

    # Calculate p-value (two-tailed)
    p_value <- mean(abs(permutation_values) >= abs(observed_value))

    # Calculate bias estimate
    bias_estimate <- observed_value

    if (verbose) {
      message("  Mean permuted value: ", mean(permutation_values))
      message("  P-value: ", p_value)
    }

    # Store results
    #results$Observed_Value[j] <- observed_value
    results$Mean_Permuted_Value[j] <- mean(permutation_values)
    results$P_Value[j] <- p_value
    results$Bias_Estimate[j] <- bias_estimate

    # Store permutation values for plotting
    permutation_values_list[[fe_name]] <- permutation_values
  }

  # Format the values in the final table
  #observed_value is bias_estimate
 # results$Observed_Value <- round(results$Observed_Value, 5)
  results$Mean_Permuted_Value <- round(results$Mean_Permuted_Value, 5)
  results$P_Value <- round(as.numeric(results$P_Value), 5)
  results$Bias_Estimate <- round(results$Bias_Estimate, 5)

  # Create result object
  result_object <- list(
    results_table = results,
    permutation_values = permutation_values_list,
    model = model
  )

  class(result_object) <- "mixedbiastest"
  return(result_object)
}
