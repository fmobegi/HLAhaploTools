# ~~~~~~~~~~~~~~~~~~~~~~~~~
# plot_hla_allele_frequency
# ~~~~~~~~~~~~~~~~~~~~~~~~~
#' Plot HLA Allele Frequencies by Gene
#'
#' Visualizes allele frequencies across genes, grouping rare alleles into "others".
#' Automatically splits classical and non-classical HLA genes into separate panels with nested facets.
#'
#' @param freq_data A data frame with gene, allele, and freq columns.
#' @param min_freq Minimum frequency threshold (default = 0.05) above which allele labels are shown.
#' @param split_threshold Threshold for number of unique genes that triggers multi-panel layout (default = 11).
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A ggplot object that can be printed or composed with patchwork.
#'
#' @import ggplot2
#' @importFrom dplyr group_by mutate summarise arrange desc filter case_when
#' @importFrom tidyr replace_na
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom ggh4x facet_nested
#'
#' @export
plot_hla_allele_frequency <- function(freq_data,
                                      min_freq = 0.05,
                                      split_threshold = 11,
                                      quiet = FALSE) {
   if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("❌ Please install ggplot2")
   }
   if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("❌ Please install patchwork")
   }
   if (!requireNamespace("ggh4x", quietly = TRUE)) {
      stop("❌ Please install ggh4x")
   }

   if (!quiet) {
      message("\t📊 Starting allele frequency plot generation...")
   }

   req_cols <- c("gene", "allele", "freq")
   if (!all(req_cols %in% colnames(freq_data))) {
      stop(
         "\t❌ Missing required columns: ",
         paste(setdiff(
            req_cols,
            colnames(freq_data)
         ), collapse = ", ")
      )
   }

   # Gene classification
   class_i <- c("A", "B", "C")
   class_ii <- c("DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5")
   nonclass_i <- c("E", "F", "G", "H", "J", "MICA", "MICB", "HFE")
   nonclass_ii <- c("DMA", "DMB", "DOA", "DOB")

   freq_data <- freq_data %>%
      dplyr::mutate(
         super_class = dplyr::case_when(
            gene %in% c(class_i, class_ii) ~ "Classical",
            gene %in% c(nonclass_i, nonclass_ii) ~ "Nonclassical",
            TRUE ~ "Other"
         ),
         labels = dplyr::case_when(
            gene %in% class_i ~ "Class I",
            gene %in% class_ii ~ "Class II",
            gene %in% nonclass_i ~ "Class I",
            gene %in% nonclass_ii ~ "Class II",
            TRUE ~ "Other"
         )
      )

   unique_genes <- unique(freq_data$gene)
   need_split <- length(unique_genes) > split_threshold

   # Show genes present by class
   classical_genes_present <- sort(
      unique(freq_data$gene[freq_data$super_class == "Classical"])
   )
   nonclassical_genes_present <- sort(
      unique(freq_data$gene[freq_data$super_class == "Nonclassical"])
   )

   if (!quiet) {
      message(
         "\t🔍 Classical genes present: ",
         paste(classical_genes_present, collapse = ", ")
      )
      message(
         "\t🔍 Nonclassical genes present: ",
         paste(nonclassical_genes_present, collapse = ", ")
      )
      message(
         sprintf(
            "\t🔢 Total unique genes found: %d (split threshold: %d)",
            length(unique_genes),
            split_threshold
         )
      )
      if (need_split) {
         message("\t⚙️ Plot will be split into Classical and Nonclassical panels.")
      } else {
         message("\t⚙️ Plotting all genes in a single panel.")
      }
   }

   # Faceted panel builder
   create_nested_panel <- function(data, title = NULL) {
      if (nrow(data) == 0) {
         ggplot2::ggplot() +
            ggplot2::theme_void()
      }

      if (!quiet) {
         message(sprintf(
            "\t\t🎨 Creating panel for %d rows of data...",
            nrow(data)
         ))
      }

      plot_data <- data %>%
         dplyr::group_by(gene) %>%
         dplyr::mutate(
            allele_label = ifelse(freq >= min_freq, allele, "Other"),
            allele_label = tidyr::replace_na(allele_label, "Other")
         ) %>%
         dplyr::group_by(super_class, labels, gene, allele_label) %>%
         dplyr::summarise(freq = sum(freq), .groups = "drop") %>%
         dplyr::arrange(gene, dplyr::desc(freq))

      gene_colors <- create_gene_palette(unique(plot_data$gene))
      color_map <- create_freq_gradient_palette(plot_data, gene_colors)

      ggplot2::ggplot(
         plot_data,
         ggplot2::aes(
            x = gene,
            y = freq,
            fill = interaction(gene, allele_label)
         )
      ) +
         ggplot2::geom_bar(
            stat = "identity",
            position = "stack",
            color = "white",
            linewidth = 0.2
         ) +
         ggplot2::geom_text(
            ggplot2::aes(label = ifelse(
               freq >= min_freq * 2,
               sprintf("%s\n(%.1f%%)", allele_label, freq * 100),
               ""
            )),
            position = ggplot2::position_stack(vjust = 0.5),
            size = 2.5
         ) +
         ggplot2::scale_fill_manual(values = color_map) +
         ggh4x::facet_nested(. ~ super_class + labels,
            scales = "free",
            space = "free"
         ) +
         ggplot2::labs(title = title, y = "Frequency", x = NULL) +
         ggplot2::theme_minimal() +
         ggplot2::theme(
            legend.position = "none",
            plot.title = ggplot2::element_text(
               hjust = 0.5,
               face = "bold",
               size = 12
            ),
            strip.text = ggplot2::element_text(face = "bold", size = 10),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
         )
   }

   if (need_split) {
      if (!quiet) {
         message("\t📑 Creating nested panels for Classical and Nonclassical genes...")
      }

      classical_data <- dplyr::filter(
         freq_data,
         super_class == "Classical"
      )
      nonclassical_data <- dplyr::filter(
         freq_data,
         super_class == "Nonclassical"
      )

      panels <- list()
      if (nrow(classical_data) > 0) {
         panels <- c(panels, list(create_nested_panel(classical_data)))
      }
      if (nrow(nonclassical_data) > 0) {
         panels <- c(panels, list(create_nested_panel(nonclassical_data)))
      }

      if (length(panels) == 2) {
         p <- panels[[1]] / panels[[2]] +
            patchwork::plot_layout(ncol = 1, heights = c(3, 2))
      } else {
         p <- panels[[1]]
      }

      p <- p + patchwork::plot_annotation(
         title = "HLA Allele Frequencies by Gene Class",
         theme = ggplot2::theme(
            plot.title = ggplot2::element_text(
               hjust = 0.5,
               size = 14,
               face = "bold"
            )
         )
      )
   } else {
      if (!quiet) {
         message("\t🖼️ Plotting all genes in one panel")
      }
      p <- create_nested_panel(freq_data, "HLA Allele Frequencies by Gene")
   }

   if (!quiet) {
      message("\t✅ Allele Frequency Plot created successfully.")
   }
   p
}
