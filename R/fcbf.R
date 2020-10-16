#' Fast Correlation-Based Filter for Feature Selection
#'
#' @author David Senhora Navega
#' @export
#'
#' @param x a symmetric numeric matrix with inter-correlation (association) of
#' the features (input) under study. See Details.
#' @param y a vector with the correlation (association) of x (features) with the
#' class (output). See Details
#' @return a logical vector of selected features
#'
#' @details
#' In this implementation, and following Hall (1999), correlation is used its
#' more general sense. Is up to the user to compute the matrix x and vector y
#' using the correlation / association measure that best suite the data and
#' problem at hand.
#'
#' @references
#' Hall, M.A., 1999. Correlation-based feature selection for machine learning.
#'
#' Yu, L. and Liu, H., 2003. Feature selection for high-dimensional data:
#' A fast correlation-based filter solution. In Proceedings of the
#' 20th international conference on machine learning (ICML-03) (pp. 856-863).
#'
fcbf <- function(x, y) {

  remove <- function(from, element) {
    from[!from %in% element]
  }

  track <- function(keep, current.order) {

    if (length(current.order) < 1)
      return(keep)

    current <- current.order[1]
    keep <- c(keep, current)
    current.order <- remove(from = current.order, element = current)
    mask.order <- sapply(current.order, function(i) {x[current, i] < y[i]})

    if (length(mask.order) == 0)
      return(keep)

    current.order <- current.order[mask.order]
    track(keep = keep, current.order = current.order)

  }

  # Algorithm
  keep <- NULL
  current.order  <- order(x = y, decreasing = T)
  keep <- track(keep = keep, current.order = current.order)

  select <- rep(x = F, time = ncol(x))
  names(select) <- colnames(x)
  select[keep] <- T

  return(select)

}
