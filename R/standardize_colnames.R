#' Standardize HLA Genotype Column Names
#'
#' Converts genotype column names (e.g. A.1) to A_1 format and applies consistent styling.
#'
#' @param df A raw HLA genotype data frame.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A cleaned data frame with standardized column names.
#' @examples
#' \dontrun{
#' raw_data <- read.csv("hla_data.csv")
#' clean_data <- standardize_colnames(raw_data)
#' }
#' @export
#' @importFrom janitor clean_names
standardize_colnames <- function(df, quiet = FALSE) {
   if (!quiet) message("\t🔧 Standardizing genotype column names...")
   if (!is.data.frame(df)) stop("\t❌ Input must be a data frame or tibble.")

   names(df) <- gsub("[\\.\\- ]", "_", names(df))
   df <- janitor::clean_names(df, case = "none")

   if (!quiet) {
      hla_cols <- sum(grepl("_[12]$", names(df)))
      message(sprintf(
         "\t✅ Standardized %d column names (%d HLA allele columns).",
         length(names(df)), hla_cols
      ))
   }
   df
}
