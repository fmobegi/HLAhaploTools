#' Plot Top HLA Haplotypes by Frequency
#'
#' Visualizes the most common haplotypes and their frequencies.
#'
#' @param hap_freq A tibble with haplotype frequencies.
#' @param n_top Number of top haplotypes to display (default = 20).
#' @param min_freq Minimum frequency to include (default = 0.01).
#' @param color Fill color for bars (default = "steelblue").
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' result <- run_haplo_em("family_typing_data.tsv", loci)
#' p <- plot_top_haplotypes(result, n_top = 10, min_freq = 0.02)
#' }
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_y_continuous theme_minimal labs element_text coord_flip theme ggtitle theme_void
#' @export
plot_top_haplotypes <- function(hap_freq,
                                n_top = 20,
                                min_freq = 0.01,
                                color = "steelblue",
                                quiet = FALSE) {
   if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("\t❌ Package ggplot2 is required. Install with install.packages('ggplot2')")
   }

   if (!quiet) message("\t📊 Creating haplotype frequency plot...")

   if (!all(c("Haplotype", "Frequency") %in% names(hap_freq))) {
      stop("\t❌ Input must have 'Haplotype' and 'Frequency' columns.")
   }

   plot_data <- hap_freq[hap_freq$Frequency >= min_freq, ]
   plot_data <- plot_data[order(plot_data$Frequency, decreasing = TRUE), ]

   if (nrow(plot_data) == 0) {
      warning(sprintf("\t⚠️ No haplotypes meet minimum frequency threshold of %.4f", min_freq))
      return(ggplot2::ggplot() +
         ggplot2::theme_void() +
         ggplot2::ggtitle("No haplotypes meet frequency threshold"))
   }

   if (nrow(plot_data) > n_top) {
      plot_data <- plot_data[seq_len(n_top), ]
      if (!quiet) {
         message(sprintf("\t✓ Showing top %d haplotypes.", n_top))
      }
   }

   plot_data$Haplotype <- factor(plot_data$Haplotype, levels = rev(plot_data$Haplotype))

   p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Haplotype, y = Frequency)) +
      ggplot2::geom_bar(stat = "identity", fill = color, width = 0.7) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", Frequency * 100)), hjust = -0.2) +
      ggplot2::scale_y_continuous(
         labels = scales::percent,
         limits = c(0, max(plot_data$Frequency) * 1.2)
      ) +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "Top HLA Haplotype Frequencies", x = "Haplotype", y = "Frequency") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
         axis.text.y = ggplot2::element_text(size = 10),
         panel.grid.major.y = ggplot2::element_blank(),
         panel.grid.minor = ggplot2::element_blank()
      )
}
