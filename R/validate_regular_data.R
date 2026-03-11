# ~~~~~~~~~~~~~~~~~~~~~
# validate_regular_data
# ~~~~~~~~~~~~~~~~~~~~~
#' Validate Regular Typing Data Structure
#'
#' Ensures that non-family typing data has proper structure before processing.
#'
#' @param df A data frame with HLA typing data
#' @param stop_if_invalid Logical. If TRUE, stops execution if data is invalid
#' @param verbose Logical. If TRUE, prints detailed validation information
#'
#' @return Logical indicating if data is valid
#' @keywords internal
validate_regular_data <- function(df,
                                  stop_if_invalid = TRUE,
                                  verbose = TRUE) {
   # Use structured detection
   detect <- detect_data_type(df, quiet = !verbose)
   is_family <- detect$is_family

   if (is_family) {
      if (stop_if_invalid) {
         stop(
            "\n\t❌ Error: Data appears to be in family study format. ",
            "Use family-specific functions for processing."
         )
      }
      return(FALSE)
   }
   if (verbose) {
      message("\n\t🔍 Validating regular typing data structure...")

      id_cols <- detect$id_cols
      id_col <- if (length(id_cols) > 0) id_cols[1] else NULL

      if (!is.null(id_col) && id_col %in% colnames(df)) {
         n_unique <- length(unique(df[[id_col]]))
         n_total <- nrow(df)
         message(sprintf(
            "\t   Found %d unique samples in %d records using column '%s'",
            n_unique, n_total, id_col
         ))

         if (n_unique < n_total) {
            message("\t⚠️ Warning: Duplicate sample IDs detected")
         }
      } else {
         message("\t⚠️ No valid sample identifier found — using row numbers")
      }
   }

   return(TRUE)
}
