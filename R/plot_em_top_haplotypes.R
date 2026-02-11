#' Plot Top 20 Haplotypes by EM Probability
#'
#' Produces a bar chart showing the top 20 haplotypes ranked by EM probability.
#' Bars are coloured by match quality, allowing visual comparison of how well
#' the highest‑probability EM haplotypes align with segregation‑based haplotypes.
#'
#' @param comparison_df A data frame containing at least the columns:
#'   \describe{
#'     \item{em_probability}{Numeric EM probability values.}
#'     \item{match_quality}{A factor describing match quality categories.}
#'   }
#'
#' @return A ggplot2 object representing the ranked haplotype probability plot.
#'
#' @importFrom dplyr arrange desc mutate n
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual labs theme_minimal
#' @importFrom ggplot2 theme element_text element_blank
#' @importFrom stats reorder
#'
#' @keywords internal
plot_em_top_haplotypes <- function(comparison_df) {
   top_haplotypes <- comparison_df %>%
      dplyr::arrange(dplyr::desc(em_probability)) %>%
      head(20) %>%
      dplyr::mutate(
         rank = 1:dplyr::n(),
         label = paste0("Rank ", rank)
      )

   ggplot(top_haplotypes, aes(
      x = reorder(label, -rank), y = em_probability,
      fill = match_quality
   )) +
      geom_col(color = "black", alpha = 0.8) +
      scale_fill_manual(
         values = c(
            "Exact Match" = "darkgreen",
            "Excellent (≥90%)" = "green3",
            "Good (75-89%)" = "orange",
            "Partial (50-74%)" = "red",
            "Poor (<50%)" = "darkred"
         )
      ) +
      labs(
         title = "Top 20 Most Common Haplotypes (by EM Probability)",
         subtitle = "Colors indicate match quality with segregation analysis",
         x = "Haplotype Rank",
         y = "EM Probability",
         fill = "Match Quality"
      ) +
      theme_minimal(base_size = 12) +
      theme(
         plot.title = element_text(face = "bold", size = 14),
         axis.text.x = element_text(angle = 45, hjust = 1),
         legend.position = "bottom",
         panel.grid.minor = element_blank()
      )
}
