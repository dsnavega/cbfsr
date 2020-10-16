#' Compute Cramer' Phi (V) between two factors
#'
#' @author David Senhora Navega
#' @import stats
#' @export
#'
#' @param  x a factor
#' @param  y a factor
#' @param  bias a logical stating if bias correction applied. Default is TRUE.
#' @param ... additional arguments to \code{\link[stats]{chisq.test}}
#'
#' @return an htest object
#'
cramer.test <- function(x, y, bias = T, ...) {

  # Pearson's Chi Square Test
  test <- suppressWarnings(expr = stats::chisq.test(x = x, y = y))

  chi <- test$statistic
  n <- sum(test$observed)
  ndim <- dim(test$observed)

  # Cramer's V from Pearson's Chi Square
  phi <- sqrt(chi / (n * min(ndim - 1)))
  names(phi) <- "V"

  if (bias) {
    # Bias Correction
    expectation <- prod(ndim - 1) / (n - 1)
    phi <- max(0, (chi / n) - expectation)
    mdim <- ndim - ((ndim - 1) ^ 2 / (n - 1))
    phi <- sqrt(phi / min(mdim - 1))
    names(phi) <- "V"
  }

  test$statistic <- phi
  test$method <- "Cramer's V"
  test$data.name <- "x (factor) & y (factor)"
  test$n <- n

  return(test)

}

#' Compute association matrix between factor variables using Cramer's Phi (V)
#'
#' @author David Senhora Navega
#' @export
#'
#' @param x a data.frame where all columns are factors.
#' @return a symmetric matrix with the pairwise phi values among the columns of
#' x
#'
phi <- function(x) {

  # Exception Handling
  if (!is.data.frame(x))
    stop("\n(-) x must be a data.frame")

  if (any(!sapply(x, is.factor)))
    stop("\n(-) All columns of x must be factors.")

  cramer.matrix <- matrix(1, nrow = ncol(x), ncol = ncol(x))
  dimnames(cramer.matrix) <- list(names(x), names(x))
  n <- ncol(x)
  for (i in seq(from = 1, to = n - 1, by = 1)) {
    for (j in seq(from = i + 1, to = n, by = 1)) {
      cramer.ij <- cramer.test(x = x[,i], y = x[,j], bias = T)$statistic
      cramer.matrix[i, j] <- round(x = cramer.ij, digits = 3)
      cramer.matrix[j, i] <- round(x = cramer.ij, digits = 3)
    }
  }

  return(cramer.matrix)

}
