#' Plot Match Rate of EM-Inferred Haplotypes
#'
#' Creates a bar plot showing the number of EM–segregation haplotype
#' comparisons that resulted in a full match versus a mismatch.
#'
#' @param comparison_df A data frame containing a logical column `match`
#'   indicating whether each EM haplotype fully matched its segregation haplotype.
#'
#' @return A ggplot2 object representing the bar plot.
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_x_discrete labs theme_minimal
#'
#' @keywords internal
plot_em_match_rate <- function(comparison_df) {
   ggplot2::ggplot(comparison_df, ggplot2::aes(x = match)) +
      ggplot2::geom_bar(fill = "orange", color = "red") +
      ggplot2::scale_x_discrete(labels = c("FALSE" = "Mismatch", "TRUE" = "Match")) +
      ggplot2::labs(
         title = "Match Rate of EM-Inferred Haplotypes",
         x = "Match Status",
         y = "Number of Haplotype Comparisons"
      ) +
      ggplot2::theme_minimal(base_size = 14)
}
