#' Clustering Correlations from Hierarchical Clustering Tree
#' @author David Senhora Navega
#' @import stats
#' @export
#'
#' @param x a matrix
#' @param y a vector
#' @param k an integer
#'
corclust <- function(x, y = NULL, k = NULL) {

  if (is.null(k))
    k <- floor(sqrt(ncol(x)))

  d <- stats::as.dist(sqrt(2 * (1  - x)))
  htree <- stats::hclust(d = d, method = "ward.D2")

  if (!is.null(y)) {
    htree <- reorder.hclust(x = htree, weights = y, agglomerator = "max")
    splitter <- cutree_order(tree = htree, k = k)
    splitter.order <- order(sapply(split(y, splitter), max), decreasing = T)
    splitter <- match(splitter, splitter.order)
    names(splitter) <- colnames(x)
  } else {
    splitter <- cutree_order(tree = htree, k = k)
  }

  object <- list(tree = htree, clust = splitter)

  return(object)

}

#' @author Jari Oksanen
#' @noRd
#' @note
#' dsnavega:This function was extracted from the vegan package.
#' Function and variables have been renamed to match my style and personal taste
#' but authorship of algorithm is given to Jari Oksanen.
#' https://github.com/vegandevs/vegan/blob/master/R/as.hclust.spantree.R
#'
hclust_mergeorder <- function(merge) {

  order <- numeric(nrow(merge) + 1)
  index <- 0
  visit <- function(i, j) {
    if (merge[i,j] < 0) {
      index <<- index + 1
      order[index] <<- -merge[i, j]
    } else {
      visit(merge[i, j], 1)
      visit(merge[i, j], 2)
    }
  }

  visit(nrow(merge), 1)
  visit(nrow(merge), 2)

  return(order)

}

#' @author Jari Oksanen
#' @noRd
#' @note
#' dsnavega:This function was extracted from the vegan package.
#' Function and variables have been renamed to match my style and personal taste
#' but authorship of algorithm is given to Jari Oksanen.
#' https://github.com/vegandevs/vegan/blob/master/R/as.hclust.spantree.R
#'
reorder.hclust <- function(x, weights,
  agglomerator = c("mean", "min", "max", "sum", "uwmean"), ...
) {

  agglomerator <- match.arg(agglomerator)
  merge <- x$merge
  nlevel <- nrow(merge)
  stats <- numeric(nlevel)
  counts <- numeric(nlevel)
  pair <- numeric(2)
  pairw <- numeric(2)

  for (i in 1:nlevel) {
    for (j in 1:2) {
      if (merge[i, j] < 0) {
        pair[j] <- weights[-merge[i, j]]
        pairw[j] <- 1
      } else {
        pair[j] <- stats[merge[i, j]]
        pairw[j] <- counts[merge[i, j]]
      }
    }
    merge[i,] <- merge[i, order(pair)]

    stats[i] <- switch(agglomerator,
      "mean" = weighted.mean(pair, pairw),
      "min" = min(pair),
      "max" = max(pair),
      "sum" = sum(pair),
      "uwmean" = mean(pair)
    )
    counts[i] <- sum(pairw)
  }

  order <- hclust_mergeorder(merge)
  x$merge <- merge
  x$order <- order
  x$value <- stats

  return(x)

}

#' @author Jari Oksanen
#' @import stats
#' @noRd
#' @note
#' dsnavega:This function was extracted from the vegan package.
#' Function and variables have been renamed to match my style and personal taste
#' but authorship of algorithm is given to Jari Oksanen.
#' https://github.com/vegandevs/vegan/blob/master/R/as.hclust.spantree.R
#'
cutree_order <- function(tree, k = NULL, h = NULL) {

  cut <- stats::cutree(tree = tree, k, h)
  if (!is.matrix(cut)) {
    cut <- order(unique(cut[tree$order]))[cut]
    names(cut) <- tree$labels
  } else {
    for (i in seq_len(ncol(cut))) {
      cut[, i] <- order(unique(cut[tree$order,i]))[cut[,i]]
    }
  }

  return(cut)

}

