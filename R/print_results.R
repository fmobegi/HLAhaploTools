# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# print.HLAhaploTools object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Print summary of HLAhaploTools result
#'
#' Custom print method for objects returned by `HLAhaploTools()`.
#' Displays a summary of key results including allele frequencies,
#' haplotype inference, family segregation, and deletion checks.
#'
#' @param x Output list returned from `HLAhaploTools()`
#' @param ... Ignored
#'
#' @export
print.HLAhaploTools <- function(x, ...) {
  cat("───────────────────────────────────────────────────\n")
  cat("                 OUTPUT SUMMARY\n")
  cat("───────────────────────────────────────────────────\n")

  # Typing data
  if (!is.null(x$typing_data)) {
    cat(sprintf("✅ Typing data:\t\t\t%d samples\n", nrow(x$typing_data)))
  } else {
    cat("⚠️  Typing data not present.\n")
  }

  # Allele frequencies
  if (!is.null(x$allele_frequencies)) {
    n_loci <- length(unique(x$allele_frequencies$gene))
    n_alleles <- length(unique(x$allele_frequencies$allele))
    cat(sprintf("✅ Allele summary:\t\t%d alleles across %d loci\n", n_alleles, n_loci))
  } else {
    cat("⚠️  Allele frequency summary not available.\n")
  }

  # Haplotype inference
  if (!is.null(x$haplotype_frequencies)) {
    cat(sprintf("✅ Haplotypes inferred:\t\t%d patterns\n", nrow(x$haplotype_frequencies)))
    cat(sprintf("   EM convergence:\t\t%s\n", ifelse(x$convergence, "Yes", "No")))
  } else {
    cat("⚠️  Haplotype inference skipped or failed.\n")
  }

  # Diplotype calls
  if (!is.null(x$top_diplotypes)) {
    cat(sprintf("✅ Diplotype assignments:\t%d subjects\n", nrow(x$top_diplotypes)))
  } else {
    cat("⚠️  Diplotype calls unavailable.\n")
  }

  # Family segregation
  if (!is.null(x$family_segregation)) {
    cat(sprintf("✅ Family segregation:\t\t%d haplotypes\n", nrow(x$family_segregation)))
  } else {
    cat("ℹ️  Family segregation not run or not applicable.\n")
  }

  # Deleted alleles
  if (!is.null(x$deleted_alleles)) {
    cat(sprintf("⚠️  Deleted alleles flagged:\t%d entries\n", nrow(x$deleted_alleles)))
  } else {
    cat("✅ No deleted alleles flagged.\n")
  }

  cat("───────────────────────────────────────────────────\n")
  cat("💡 Use `names(x)` to inspect available components.\n")
  cat("───────────────────────────────────────────────────\n")
  invisible(x)
}
