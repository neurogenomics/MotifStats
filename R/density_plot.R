#' Plot the positional distribution of motifs relative to peak summits
#'
#' \code{density_plot} creates a density plot of the positional distribution
#' of motifs within peaks.
#'
#' @import ggplot2
#'
#' @returns NULL
#' @export
density_plot <- function(distance_vec,
                         plot_title = NULL,
                         x_label = NULL,
                         y_label = NULL,
                         limits = c(-400, 400)){
  df <- as.data.frame(distance_vec)
  colnames(df) <- "dist"

  ggplot(df, aes(x = dist)) +
    geom_line(stat = "density", linetype = "solid", linewidth = 1) +
    labs(x = x_label,
         y = y_label,
         title = plot_title) +
    theme_minimal() +
    theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
          plot.title = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 14, margin = margin(t = 10), colour = "black"),
          axis.title.y = element_text(size = 14, margin = margin(r = 10), colour = "black"),
          axis.text = element_text(size = 11, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(color = "black", linewidth = 0.5)) +
    expand_limits(x = 0, y = -1e-04) +
    scale_x_continuous(expand = c(0, 0),
                       limits = limits,
                       breaks = seq(-1e+05, 1e+05, by = 100)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
}
