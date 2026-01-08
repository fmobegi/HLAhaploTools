# ~~~~~~~~~~~~~~~~~~~
# reorder_drb_alleles
# ~~~~~~~~~~~~~~~~~~~
#' Reorder DRB Alleles to Ensure _1 Is Populated
#'
#' Internal helper function that checks the DRB1/3/4/5 loci and ensures that
#' when an allele is only present in the `_2` column (and `_1` is missing), it is
#' moved to `_1` and `_2` is cleared. This ensures consistency and simplifies
#' downstream processing and compliance checks.
#'
#' @param df A data frame containing DRB1/3/4/5 allele columns (e.g., DRB1_1, DRB1_2).
#' @param quiet Logical. If `TRUE`, suppresses status messages.
#'
#' @return A modified data frame with reordered DRB allele columns where applicable.
#' @keywords internal
reorder_drb_alleles <- function(df, quiet = FALSE) {
  drb_genes <- c("DRB1", "DRB3", "DRB4", "DRB5")
  swapped_total <- 0 # Counter for how many swaps were made

  for (gene in drb_genes) {
    a1 <- paste0(gene, "_1")
    a2 <- paste0(gene, "_2")

    # Proceed only if both _1 and _2 columns exist
    if (all(c(a1, a2) %in% names(df))) {
      for (i in seq_len(nrow(df))) {
        allele1 <- df[[a1]][i]
        allele2 <- df[[a2]][i]

        # Swap if _1 is empty/missing and _2 is populated
        if ((is.na(allele1) || allele1 == "") &&
          !(is.na(allele2) || allele2 == "")) {
          df[[a1]][i] <- allele2
          df[[a2]][i] <- NA
          swapped_total <- swapped_total + 1
        }
      }
    }
  }

  if (!quiet && swapped_total > 0) {
    message(sprintf(
      "\t🔄 Reordered DRB alleles in %d samples to ensure _1 is populated when possible.",
      swapped_total
    ))
  }
  df
}
