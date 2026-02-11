#' Plot EM vs Segregation Haplotype Match Distribution
#'
#' Generates a histogram showing the distribution of percentage match
#' between EM and segregation haplotypes.
#'
#' @param comparison_df A data frame containing a column named
#'   `percentage_match` with numeric match percentages.
#'
#' @return A ggplot2 object representing the histogram.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram scale_x_continuous labs theme_minimal
#'
#' @keywords internal
#'
plot_em_segregation_distribution <- function(comparison_df) {
   ggplot2::ggplot(comparison_df, ggplot2::aes(x = percentage_match)) +
      ggplot2::geom_histogram(
         binwidth = 5,
         boundary = 0,
         closed = "left",
         fill = "orange",
         color = "red"
      ) +
      ggplot2::scale_x_continuous(
         limits = c(40, 100),
         breaks = seq(40, 100, by = 10)
      ) +
      ggplot2::labs(
         title = "Distribution of EM vs Segregation Haplotype Match",
         x = "Percentage Match (%)",
         y = "Number of EM–Segregation Matches"
      ) +
      ggplot2::theme_minimal(base_size = 14)
}
