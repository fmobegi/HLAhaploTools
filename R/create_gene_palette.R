#' Create a palette of distinct colors for HLA genes
#'
#' Generates a named vector of colors for a set of HLA genes.
#' Known genes receive predefined colors; unknown genes are assigned
#' dynamically generated colors using golden‑ratio HSV spacing.
#'
#' @param genes Character vector of gene names.
#'
#' @return A named character vector of hex colors.
#'
#' @keywords internal
#' @noRd
create_gene_palette <- function(genes) {
  hla_colors <- c(
    "A" = "#E41A1C", "B" = "#377EB8", "C" = "#4DAF4A",
    "DRB1" = "#984EA3", "DRB3" = "#925E93", "DRB4" = "#7F4E83",
    "DRB5" = "#6C3E73", "DQB1" = "#FF7F00", "DQA1" = "#FFFF33",
    "DPB1" = "#A65628", "DPA1" = "#F781BF", "E" = "#999999",
    "F" = "#66C2A5", "G" = "#FC8D62", "H" = "#8DA0CB",
    "J" = "#E78AC3", "MICA" = "#A6D854", "MICB" = "#FFD92F",
    "HFE" = "#B3B3B3", "DMA" = "#B2DF8A", "DMB" = "#FB9A99",
    "DOA" = "#CAB2D6", "DOB" = "#FDBF6F"
  )

  result <- character(length(genes))
  names(result) <- genes
  dynamic_index <- 1

  for (i in seq_along(genes)) {
    gene <- genes[i]
    if (gene %in% names(hla_colors)) {
      result[i] <- hla_colors[gene]
    } else {
      h <- (dynamic_index * 0.618033988749895) %% 1
      result[i] <- hsv(h, 0.6, 0.9)
      dynamic_index <- dynamic_index + 1
    }
  }

  result
}
