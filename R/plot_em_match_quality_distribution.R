#' Plot EM Match Quality Distribution
#'
#' Creates a bar chart summarising the distribution of match quality categories
#' for EM–segregation haplotype comparisons. The function expects a data frame
#' containing a `match_quality` column and computes counts and percentages
#' before generating the plot.
#'
#' @param comparison_df A data frame containing a factor column `match_quality`
#'   with levels such as "Exact Match", "Excellent (≥90%)", "Good (75-89%)",
#'   "Partial (50-74%)", and "Poor (<50%)".
#'
#' @return A ggplot2 object representing the match quality distribution plot.
#'
#' @importFrom dplyr group_by summarise mutate n
#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_fill_manual labs
#' @importFrom ggplot2 theme_minimal theme element_text element_blank ylim
#'
#' @keywords internal
plot_em_match_quality_distribution <- function(comparison_df) {
   # Defensive check: ensure match_quality exists
   if (!"match_quality" %in% names(comparison_df)) {
      stop("comparison_df must contain a 'match_quality' column.")
   }

   match_summary <- comparison_df %>%
      dplyr::group_by(match_quality) %>%
      dplyr::summarise(
         count = dplyr::n(),
         .groups = "drop"
      ) %>%
      dplyr::mutate(
         percentage = (count / sum(count)) * 100,
         label = paste0(count, "\n(", round(percentage, 1), "%)")
      )

   ggplot(match_summary, aes(x = match_quality, y = count, fill = match_quality)) +
      geom_col(color = "black", alpha = 0.8) +
      geom_text(aes(label = label), vjust = -0.5, size = 4, fontface = "bold") +
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
         title = "Match Quality Distribution",
         subtitle = "EM haplotypes categorized by match percentage with segregation",
         x = "Match Quality Category",
         y = "Number of Haplotypes"
      ) +
      theme_minimal(base_size = 12) +
      theme(
         plot.title = element_text(face = "bold", size = 14),
         legend.position = "none",
         axis.text.x = element_text(angle = 0, hjust = 0.5),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank()
      ) +
      ylim(0, max(match_summary$count) * 1.15)
}
