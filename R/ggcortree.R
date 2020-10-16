#' Visualize a Correlation Matrix with an Hierarchical Tree using ggplot2
#'
#' @author David Senhora Navega
#' @export
#'
#' @import stats ggdendro ggplot2
#'
#' @param x a matrix of pairwise correlation among features
#' @param y vector of correlation of the features represented in x with a target
#' quantity. Default is NULL which means that this information is not exploited.
#' If y is supplied it is used to reorder the hierarchical tree under the
#' constrains of its structure.
#' @param k an integer. Number of cluster to cut the tree into. Default is NULL.
#' k is estimate as the minimal number of clusters the would result in an isolated
#' feature.
#' @param circular should the tree be represent in polar coordinates? Default F
#' @param nudge a numeric.  Displacement factor of the labels. Aesthetics purpose.
#' @param expand a numeric. Expansion of the y-axis. Aesthetics purpose.
#' @param label.size a numeric. Size of the terminal nodes of the tree labels.
#'
ggcortree <- function(x,
  y = NULL, k = NULL, circular = F, nudge = 0.05, expand = 0.2, label.size = 3.5
) {

  n <- nrow(x)

  x <- (1 - x)
  htree <- stats::hclust(d = as.dist(x), method = "complete")

  if (!is.null(y))
    htree <- reorder.hclust(x = htree, weights = y ^ 2, agglomerator = "max")

  if (is.null(k)) {

    kclust.test <- sapply(2:(n - 1), function(i) {
      test <- min(table(cutree_order(tree = htree, k = i))) > 1
      return(test)
    })
    k <- sum(kclust.test) + 1

  }


  hcdata <- dendro_data_k(hc = htree, k = k)
  if (is.null(y)) {
    hcdata$labels$size <- rep(x = 1, times = n)
  } else {
    y <- y ^ 2
    hcdata$labels$size <- y[match(hcdata$labels$label, names(y))]
  }

  size <- NULL
  p <- ggdendrogram(
    hcdata = hcdata,
    direction = "rl",
    circular = circular,
    branch.size = 1,
    label.size = label.size,
    expand.y = expand,
    nudge.label = nudge
  ) + ggplot2::geom_point(
    data = hcdata$labels,
    mapping = ggplot2::aes(
      x = match(label, label),
      y = 0,
      size = size * 0.5
    ),
    shape = 16,
    show.legend = F
  )

  return(p)

}

#' @author atrebas.github.io
#' @noRd
#' @importFrom ggdendro segment label
#' @note
#' dsnavega: This function and/or variables might been renamed to match
#' my personal code style but key functionality was implement by github user
#' atrebas
#' https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
#'
ggdendrogram <- function(hcdata,
  circular = FALSE,
  direction = c("lr", "rl", "tb", "bt"),
  branch.size = 1,
  label.size  = 3,
  nudge.label = 0.05,
  expand.y = 0.1
) {

  direction <- match.arg(direction)
  ybreaks <- pretty(ggdendro::segment(hcdata)$y, n = 5)
  ymax <- max(ggdendro::segment(hcdata)$y)

  # Branches
  x <- y <- xend <- yend <- angle <- size <- NULL
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = ggdendro::segment(hcdata),
      mapping = ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        linetype = factor(line)
      ),
      lineend = "round", show.legend = FALSE, size = branch.size
    )

  # Orientation
  if (circular) {
    p <- p +
      ggplot2::coord_polar(direction = -1) +
      ggplot2::scale_x_continuous(
        breaks = NULL, limits = c(0, nrow(ggdendro::label(hcdata)))
      ) +
      ggplot2::scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p +
      ggplot2::scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + ggplot2::coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p +
        ggplot2::scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p +
        ggplot2::scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }

  # Labels
  label.params <- set_labels_params(nrow(hcdata$labels), direction, circular)
  hcdata$labels$angle <- label.params$angle

  p <- p +
    ggplot2::geom_text(data = ggdendro::label(hcdata),
      mapping = ggplot2::aes(
        x = x,
        y = y,
        label = label,
        angle = angle
      ),
      vjust = label.params$vjust,
      hjust = label.params$hjust,
      nudge_y = ymax * nudge.label,
      size = label.size,
      show.legend = FALSE
    )

  # Limits
  ylim <- -round(ymax * expand.y, 1)
  p <- p +
    ggplot2::expand_limits(y = ylim)

  # Theme
  p <- p +
    ggplot2::theme_void()

  return(p)

}

#' @author atrebas.github.io
#' @noRd
#'
#' @importFrom ggdendro dendro_data
#' @import stats
#' @note
#' dsnavega: This function name or variables might been renamed to match
#' my personal code style but key functionality was implement by github user
#' atrebas
#' https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
#'
dendro_data_k <- function(hc, k) {

  hcdata <- ggdendro::dendro_data(hc, type = "rectangle")
  seg <- hcdata$segments
  labclust <- stats::cutree(hc, k)[hc$order]
  segclust <- rep(0L, nrow(seg))
  heights <- sort(hc$height, decreasing = TRUE)
  height <- mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)

  for (i in 1:k) {
    xi <- hcdata$labels$x[labclust == i]
    idx1 <- seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2 <- seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3 <- seg$yend < height
    idx <- idx1 & idx2 & idx3
    segclust[idx] <- i
  }

  idx <- which(segclust == 0L)
  segclust[idx] <- segclust[idx + 1L]
  hcdata$segments$clust <- segclust
  hcdata$segments$line <- as.integer(segclust < 1L)
  hcdata$labels$clust <- labclust

  return(hcdata)

}

#' @author atrebas.github.io
#' @noRd
#' @note
#' dsnavega: This function name or variables might been renamed to match
#' my personal code style but key functionality was implement by github user
#' atrebas.
#' https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
#'
set_labels_params <- function(
  nblabels, direction = c("tb", "bt", "lr", "rl"), circular = FALSE
) {

  if (circular) {
    angle <- 360 / nblabels * 1:nblabels + 90
    idx <- angle >= 90 & angle <= 270
    angle[idx] <- angle[idx] + 180
    hjust <- rep(0, nblabels)
    hjust[idx] <- 1
  } else {
    angle <- rep(0, nblabels)
    hjust <- 0
    if (direction %in% c("tb", "bt")) {angle <- angle + 45}
    if (direction %in% c("tb", "rl")) {hjust <- 1}
  }

  params <- list(angle = angle, hjust = hjust, vjust = 0.5)
  return(params)

}

#' @author Jari Oksanen
#' @noRd
#' @note
#' dsnavega:This functions was extracted from the vegan package.
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
#' dsnavega:This functions was extracted from the vegan package.
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
#' dsnavega:This functions was extracted from the vegan package.
#' Function and variables have been renamed to match my style and personal taste
#' but authorship of algorithm is given to Jari Oksanen.
#' https://github.com/vegandevs/vegan/blob/master/R/as.hclust.spantree.R
#'
cutree_order <- function(tree, k = NULL, h = NULL) {

  cut <- stats::cutree(tree, k, h)
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
