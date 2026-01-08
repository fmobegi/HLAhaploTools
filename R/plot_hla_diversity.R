# ~~~~~~~~~~~~~~~~~~
# plot_hla_diversity
# ~~~~~~~~~~~~~~~~~~
#' Plot HLA Allele Diversity by Population
#'
#' Shows the relative frequency of top alleles for a specific gene within each population.
#'
#' @param freq_data A data frame with gene, allele, freq and pop columns.
#' @param gene The HLA gene to analyze (e.g., "A", "DRB1").
#' @param ntop Number of top alleles to label per population (default = 5).
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A ggplot object that can be further customized or printed.
#'
#' @import ggplot2
#' @importFrom forcats fct_lump_n
#' @importFrom dplyr group_by mutate summarise arrange desc
#'
#' @export
plot_hla_diversity <- function(freq_data,
                               gene = "A",
                               ntop = 5,
                               quiet = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "\t❌ Package 'ggplot2' is required. Please install it with install.packages('ggplot2')"
    )
  }
  if (!requireNamespace("forcats", quietly = TRUE)) {
    stop(
      "\t❌ Package 'forcats' is required. Please install it with install.packages('forcats')"
    )
  }

  if (!quiet) {
    message(sprintf("\t📊 Generating allele diversity plot for gene %s...", gene))
  }

  # Validate input data
  req_cols <- c("gene", "allele", "freq", "pop")
  if (!all(req_cols %in% colnames(freq_data))) {
    stop("\t❌ Input data must contain columns: gene, allele, freq, pop")
  }

  # Filter data for the selected gene
  gene_data <- freq_data[freq_data$gene == gene, ]

  if (nrow(gene_data) == 0) {
    warning(sprintf("\t⚠️ No allele data found for gene '%s'.", gene))

    ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::ggtitle(sprintf("No allele data for %s", gene))
  }

  if (!quiet) {
    message("\t🔍 Filtering top ", ntop, " alleles per population...")
  }

  # For each population, keep only top N alleles and group the rest as "Other"
  plot_data <- gene_data %>%
    dplyr::group_by(pop) %>%
    dplyr::mutate(allele_group = forcats::fct_lump_n(
      allele,
      n = ntop,
      w = freq,
      other_level = "Other"
    )) %>%
    dplyr::group_by(pop, allele_group) %>%
    dplyr::summarise(freq = sum(freq), .groups = "drop") %>%
    dplyr::arrange(pop, dplyr::desc(freq))

  # Generate color palette specific to population+allele combinations
  color_map <- create_pop_allele_palette(plot_data, gene)

  # Compose the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = pop,
    y = freq,
    fill = interaction(pop, allele_group)
  )) +
    ggplot2::geom_bar(
      stat = "identity",
      position = "stack",
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(
        freq >= 0.05,
        sprintf("%s\n(%.1f%%)", allele_group, freq * 100),
        ""
      )),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 2.8
    ) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::labs(
      title = sprintf("%s Allele Distribution by Population", gene),
      x = "Population",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(
        face = "bold",
        angle = 45,
        hjust = 1
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 13,
        face = "bold"
      )
    )

  if (!quiet) {
    message("\t✅ Diversity plot created successfully.")
  }
  p
}
