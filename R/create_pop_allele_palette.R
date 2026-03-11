#' Create population‑specific allele gradient palettes
#'
#' Generates a gradient palette for each allele within each population,
#' based on a single gene's base color. Lower‑frequency alleles fade
#' more strongly toward white.
#'
#' @param data A data frame containing `pop`, `freq`, and `allele_group`.
#' @param gene Character string giving the gene name.
#'
#' @return A named character vector of hex colors for each pop–allele ID.
#'
#' @keywords internal
#' @noRd
create_pop_allele_palette <- function(data, gene) {
   gene_colors <- create_gene_palette(gene)
   base_color <- gene_colors[[1]]
   result <- character()

   for (pop in unique(data$pop)) {
      pop_data <- data[data$pop == pop, ]
      pop_data <- pop_data[order(pop_data$freq, decreasing = TRUE), ]
      n_alleles <- nrow(pop_data)
      rgb_base <- col2rgb(base_color) / 255

      for (i in 1:n_alleles) {
         fade <- min(0.8, (i - 1) / max(1, n_alleles - 1))
         r <- rgb_base[1, 1] + (1 - rgb_base[1, 1]) * fade
         g <- rgb_base[2, 1] + (1 - rgb_base[2, 1]) * fade
         b <- rgb_base[3, 1] + (1 - rgb_base[3, 1]) * fade
         color <- rgb(r, g, b)

         id <- paste(pop, pop_data$allele_group[i], sep = ".")
         result[id] <- color
      }
   }

   result
}
