#' Fast Computation of Pearson and Spearman Correlation Coefficient
#'
#' @author David Senhora Navega
#' @export
#'
#' @param x a matrix or data.frame. Data is internally coerced to numeric
#' representation.
#' @param y vector (numeric or factor). Internally coerced to numeric
#' representation.
#' @param method a character stating the method to used. Default is "pearson".
#' @param partial a logical stating if correlation among x should be account
#' for y. Default is T.
#' @param weights a numeric vector of weights used to computed weighted
#' correlations.
#' @param digits an integer defining the precision of round function. Default is
#' 3
#'
#' @return a list with two components, x and y, storing the inter-correlations
#' among x and the correlations of x with y respectively. The component x
#' is a symmetric squared matrix and the component y is a numeric vector.
#'
#' @details
#' This implementation uses matrix operations for fast computation of correlation
#' coefficients. A drawback of this approach is that the data must be complete,
#' that is no NA values. This function implements a fast mean imputation step
#' where the NA values are substituted by 0 after the data been scaled as part of
#' the algorithm for fast computation of the correlation.
#'
fcbw <- function(x, y,
  method = c("pearson", "spearman"),
  partial = T,
  weights = NULL,
  digits = 3
) {

  # Parameters
  n <- nrow(x)
  m <- ncol(x)
  data <- data.frame(x, Y = y)
  xnames <- colnames(x)

  # Fast Correlation Computation as Matrix Operations
  data <- data.matrix(data)

  method <- match.arg(arg = method, choices = c("pearson", "spearman"))

  if (method == "spearman")
    data <- apply(data, 2, rank)

  if (is.null(weights))
    weights <- rep(x = 1, times = n)

  weights <- weights / sum(weights)

  data <- scale(
    x = data,
    center = apply(data, 2, weighted_mean, weights = weights),
    scale = apply(data, 2, weighted_sd, weights = weights)
  )

  data[is.na(data)] <- 0
  R <- (t(data) %*% diag(weights) %*% data)

  X <- R[-(m + 1), -(m + 1)]
  Y <- R[, m + 1][ -(m + 1)]

  # Compute Partial Correlation Using Recursive Formula as Matrix Operations
  if (partial) {

    XXY <- outer(
      X = R[, m + 1][-(m + 1)],
      Y = R[, m + 1][-(m + 1)],
      FUN = "*"
    )

    XX <- outer(
      X = sqrt(1 - R[, m + 1][-(m + 1)] ^ 2),
      Y = sqrt(1 - R[, m + 1][-(m + 1)] ^ 2),
      FUN = "*"
    )

    X <- (X - XXY) / XX

  }

  dimnames(X) <- list(xnames, xnames)
  names(Y) <- xnames

  object <- list(
    x = round(x = X, digits  = digits),
    y = round(x = Y, digits  = digits)
  )

  return(object)

}

#' @author David Senhora Navega
#' @noRd
#'
weighted_mean <- function(x, weights, na.rm = T) {
  value <- sum(x * weights, na.rm = na.rm) / sum(weights, na.rm = na.rm)
  return(value)
}

#' @author David Senhora Navega
#' @noRd
#'
weighted_sd <- function(x, weights, na.rm = T) {
  ss <- ((x - weighted_mean(x = x,weights = weights)) ^ 2)
  wss <- sum(x = ss * weights, na.rm = na.rm) / sum(x = weights, na.rm = na.rm)
  value <- sqrt(wss)
  return(value)
}
