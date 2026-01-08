#' Create gradient color palettes for gene–allele combinations
#'
#' Produces a gradient of colors for each allele of a gene, fading from
#' the gene's base color toward white according to allele frequency.
#'
#' @param data A data frame containing at least `gene`, `freq`, and
#'   either `allele_label` or `allele_group`.
#' @param gene_colors Named vector of base colors for each gene.
#'
#' @return A named character vector of hex colors for each gene–allele ID.
#'
#' @keywords internal
#' @noRd
create_freq_gradient_palette <- function(data, gene_colors) {
  result <- character()

  for (gene in unique(data$gene)) {
    gene_data <- data[data$gene == gene, ]
    base_color <- gene_colors[gene]
    gene_data <- gene_data[order(gene_data$freq, decreasing = TRUE), ]
    n_alleles <- nrow(gene_data)
    rgb_base <- col2rgb(base_color) / 255

    for (i in 1:n_alleles) {
      fade <- min(0.8, (i - 1) / max(1, n_alleles - 1))
      r <- rgb_base[1, 1] + (1 - rgb_base[1, 1]) * fade
      g <- rgb_base[2, 1] + (1 - rgb_base[2, 1]) * fade
      b <- rgb_base[3, 1] + (1 - rgb_base[3, 1]) * fade
      color <- rgb(r, g, b)

      id <- if ("allele_label" %in% colnames(gene_data)) {
        paste(gene, gene_data$allele_label[i], sep = ".")
      } else if ("allele_group" %in% colnames(gene_data)) {
        paste(gene_data$pop[i], gene_data$allele_group[i], sep = ".")
      } else {
        gene
      }

      result[id] <- color
    }
  }

  result
}
