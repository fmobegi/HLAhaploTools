#' Extract HLA Locus Names from Genotype Table
#'
#' Identifies unique locus names based on "_1" and "_2" suffixes.
#'
#' @param df A data frame with standardized column names.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A character vector of distinct locus names (e.g. A, B, DRB1).
#' @examples
#' \dontrun{
#' loci <- extract_loci(hla_data)
#' print(loci)
#' }
#' @export
extract_loci <- function(df, quiet = FALSE) {
   if (!quiet) message("\t🔬 Extracting HLA loci from column names...")
   if (!is.data.frame(df)) stop("\t❌ Input must be a data frame or tibble.")

   allele_cols <- grep("_[12]$", names(df), value = TRUE)
   if (length(allele_cols) == 0) {
      warning("\t⚠️ No HLA allele columns found with _1 or _2 suffix.")
      return(character(0))
   }

   loci <- unique(gsub("_[12]$", "", allele_cols))
   if (!quiet) {
      message(sprintf(
         "\t✅ Found %d HLA loci: %s",
         length(loci), paste(loci, collapse = ", ")
      ))
   }
   loci
}
