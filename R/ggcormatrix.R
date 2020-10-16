#' Visualize a Correlation Matrix using ggplot2
#'
#' @author David Senhora Navega
#' @export
#'
#' @import reshape2 stats ggplot2
#'
#' @param x a pairwise correlation matrix
#' @param determination logical stating if correlations should be transformed to
#' coefficients of determination by squaring the values of x.
#' @param hc.order a logical stating if rows and columns of x should be ordered
#' according the a hierarchical tree with complete linkage.
#' @param fill a user defined character vector of three strings of hexcodes
#' defining the color scheme of the legend of the plot. Default NULL. See Details.
#' @param midpoint the midpoint value for the transition of the fill gradient.
#' @param limit limit values of the legend of the plot. Default NULL. See Details.
#' @param breaks numeric values defining the ticks of the legend scale.
#' @param text.size a numeric defining the size of the plot text.
#'
#' @details
#' If determination is TRUE the color scheme of the plot is defined as:
#' fill = c(low = "#FFFFFF", mid = "#808080", high = "#232323"). Otherwise,
#' is defined as:
#' fill <- c(low = "#3B9AB2", mid = "#FFFFFF", high = "#F21A00").
#' The first option generates a gray-scale plot. The second option generates
#' a color gradient that goes from a blueish (negative correlation) to reddish
#' (positive correlation) color. In both cases a zero correlation or determination
#' coefficient is represent with white color.
#' By default the limit argument is NULL which means that this values is determined
#' internally. limit is definied a numeric vector of two elements defining the range
#' of the legend scale. For correlation it is c(-1, 1) and for determination = T
#' is c(0, 1). User can supply its own values which might be help for if observed
#' values of x have a limited or very similar range.
#'
ggcormatrix <- function(x,
  determination = T,
  hc.order = F,
  fill = NULL,
  limit = NULL,
  midpoint = NULL,
  breaks = c(0, 0.2, 0.5, 0.8, 1),
  text.size = 9
) {

  if (hc.order) {
    hc.tree <- stats::hclust(as.dist((1 - x) / 2), method = "complete")
    x <- x[hc.tree$order, hc.tree$order]
  }

  if (determination)
    x <- x ^ 2

  tile <- reshape2::melt(x, na.rm = TRUE)
  colnames(tile) <- c("x", "y", "z")

  breaks <- unique(sort(c(-breaks, breaks)))

  if (determination)
    breaks <- unique(sort(unique(breaks) ^ 2))[-2]

  if (is.null(limit))
    limit <- range(breaks)

  if (is.null(midpoint))
    midpoint <- sum(limit) / 2


  if (is.null(fill))
    if (determination) {
      fill <- c(low = "#FFFFFF", mid = "#808080", high = "#232323")
    } else {
      fill <- c(low = "#3B9AB2", mid = "#FFFFFF", high = "#F21A00")
    }

  if (determination) {
    guide.title <- "Coefficient of Determination"
  } else {
    guide.title <- "Coefficient of Correlation"
  }

  p <- ggplot2::ggplot(data = tile) +
    ggplot2::aes_string(x = "x", y = "y", fill = "z") +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = fill[1], mid = fill[2], high = fill[3],
      midpoint = midpoint, limit = limit, breaks = breaks, space = "Lab",
      guide = ggplot2::guide_colorbar(
        direction = "horizontal", title.position = "top", nbin = 100,
        ticks.colour = "white", ticks.linewidth = 2, title.hjust = 0.5,
        title = guide.title, barwidth = 16, barheight = 1
      )
    ) +
    ggplot2::theme(
      legend.position = "top",
      panel.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45, vjust = 1, hjust = 1, size = text.size - 2,
      ),
      axis.text.y = ggplot2::element_text(size = text.size)
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "", y = "")
  return(p)

}
