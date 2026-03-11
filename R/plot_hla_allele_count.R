# ~~~~~~~~~~~~~~~~~~~~~
# plot_hla_allele_count
# ~~~~~~~~~~~~~~~~~~~~~
#' Plot Count of Unique Alleles per HLA Gene
#'
#' Visualizes the number of unique observed alleles for each gene.
#'
#' @param freq_data A data frame with gene and allele columns.
#' @param min_count Minimum unique allele count to include gene (default = 0).
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A ggplot object that can be further customized or printed.
#'
#' @import ggplot2
#' @import dplyr

#' @export
plot_hla_allele_count <- function(freq_data,
                                  min_count = 0,
                                  quiet = FALSE) {
   if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop(
         "\t❌ Package ggplot2 is required for plotting. Please install it with install.packages('ggplot2')"
      )
   }

   if (!quiet) {
      message("\t📊 Starting allele count plot generation...")
   }

   # Validate input data
   req_cols <- c("gene", "allele")
   if (!all(req_cols %in% colnames(freq_data))) {
      stop("\t❌ Input data must contain columns: gene, allele")
   }

   # Gene classification to report classes found
   class_i <- c("A", "B", "C")
   class_ii <- c("DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5")
   nonclass_i <- c("E", "F", "G", "H", "J", "MICA", "MICB", "HFE")
   nonclass_ii <- c("DMA", "DMB", "DOA", "DOB")

   # Count unique alleles per gene
   if (!quiet) {
      message("\t🔍 Counting unique alleles per gene...")
   }
   allele_counts <- freq_data %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(
         unique_alleles = dplyr::n_distinct(allele),
         .groups = "drop"
      ) %>%
      dplyr::filter(unique_alleles >= min_count) %>%
      dplyr::arrange(dplyr::desc(unique_alleles))

   if (nrow(allele_counts) == 0) {
      warning(sprintf("\t⚠️ No genes with at least %d unique alleles.", min_count))
      ggplot2::ggplot() +
         ggplot2::theme_void() +
         ggplot2::ggtitle(sprintf(
            "No genes with ≥%d unique alleles", min_count
         ))
   }

   # Determine classes found in data
   found_genes <- allele_counts$gene
   classical_genes <- sort(unique(
      found_genes[found_genes %in%
         c(class_i, class_ii)]
   ))
   nonclassical_genes <- sort(unique(
      found_genes[found_genes %in%
         c(nonclass_i, nonclass_ii)]
   ))
   other_genes <- sort(unique(
      found_genes[!found_genes %in%
         c(class_i, class_ii, nonclass_i, nonclass_ii)]
   ))

   if (!quiet) {
      message("\t🧬 Genes found in each class:")
      message("\t\tClassical genes: ", ifelse(
         length(classical_genes) > 0,
         paste(classical_genes, collapse = ", "),
         "None"
      ))
      message(
         "\t\tNonclassical genes: ",
         ifelse(
            length(nonclassical_genes) > 0,
            paste(nonclassical_genes, collapse = ", "),
            "None"
         )
      )
      message("\t\tOther genes: ", ifelse(
         length(other_genes) > 0,
         paste(other_genes, collapse = ", "),
         "None"
      ))
   }

   if (!quiet) {
      message(sprintf(
         "\t🎨 Preparing plot for %d genes...",
         nrow(allele_counts)
      ))
   }
   gene_colors <- create_gene_palette(allele_counts$gene)

   p <- ggplot2::ggplot(
      allele_counts,
      ggplot2::aes(
         x = reorder(gene, unique_alleles),
         y = unique_alleles,
         fill = gene
      )
   ) +
      ggplot2::geom_bar(
         stat = "identity",
         color = "white",
         linewidth = 0.2
      ) +
      ggplot2::geom_text(
         ggplot2::aes(label = unique_alleles),
         position = ggplot2::position_stack(vjust = 0.5),
         color = "white"
      ) +
      ggplot2::scale_fill_manual(values = gene_colors) +
      ggplot2::ylab("Number of Unique Alleles") +
      ggplot2::xlab("HLA Gene") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
         legend.position = "none",
         axis.text.x = ggplot2::element_text(
            face = "bold",
            angle = 45,
            hjust = 1
         ),
         panel.grid.major = ggplot2::element_blank(),
         panel.grid.minor = ggplot2::element_blank(),
         plot.title = ggplot2::element_text(hjust = 0.5)
      ) +
      ggplot2::ggtitle("HLA Allelic Diversity")

   if (!quiet) {
      message("\t✅ Allele count plot created successfully.")
   }
   p
}
