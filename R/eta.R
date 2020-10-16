#' Compute Pearson's Eta Squared between a factor (x) and a numeric (y) vector.
#'
#' @author David Senhora Navega
#' @import stats
#' @export
#'
#' @param  x a factor vector
#' @param  y a numeric vector
#'
#' @return an htest object
#'
etasq.test <- function(x, y) {

  # Exception Handling
  if (!is.factor(x))
    stop("\n(-) x must be a factor")

  if (!(is.vector(y) & is.numeric(y)))
    stop("\n(-) x must be a factor")

  df <- na.omit(data.frame(x = x, y = y))
  for (name in c("x","y")) {
    assign(name, df[[name]])
  }

  k <- nlevels(x)
  n <- nrow(df)
  n_x <- table(x)
  y_x <- split(y, x)

  numerator <- vector(length = k)
  denominator <- vector(length = k)
  for (i in seq_len(k)) {
    numerator[i] <- (n_x[i]) * (mean(y_x[[i]]) - mean(y)) ^ 2
    denominator[i] <- sum((y_x[[i]] - mean(y)) ^ 2)
  }

  estimate <- sum(numerator) / sum(denominator)
  statistic <- n * estimate
  p.value <- pchisq(q = statistic, df = k - 1, lower.tail = F)

  frequency <- table(x)
  names(frequency) <- levels(x)
  detail <- numerator / denominator
  names(detail) <- levels(x)

  names(estimate) <- "Eta Squared"
  names(statistic) <- "Chi Squared"

  test <- structure(
    .Data = list(
      parameter = estimate, statistic = statistic, p.value = p.value,
      detail = detail, frequency = frequency, n = n, k = k, df = k - 1,
      method = "Pearson's Eta Squared",
      data.name = "x (factor) and y (numeric)"
    ),
    class = c("htest", "etasq.test")
  )

  return(test)

}


#' Compute Pearson's Eta between a data.frame of factors and a numeric vector
#'
#' @author David Senhora Navega
#' @export
#'
#' @param x a data.frame where all columns are factors.
#' @param y a numeric vector.
#' @return a vector of eta values
#'
eta <- function(x, y) {

  # Exception Handling
  if (!is.data.frame(x))
    stop("\n(-) x must be a data.frame")

  if (any(!sapply(x, is.factor)))
    stop("\n(-) All columns of x must be factors.")

  value <- sapply(x, function(x) {
    unname(sqrt(etasq.test(x = x, y = y)$parameter))
  })

  return(value)

}
