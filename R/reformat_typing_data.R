# ~~~~~~~~~~~~~~~~~~~~
# reformat_typing_data
# ~~~~~~~~~~~~~~~~~~~~
#' Reformat and Verbosely Validate Raw HLA Typing Data
#'
#' Cleans and formats HLA typing data by standardizing allele notations,
#' filling missing second allele calls with DRB1–DRB3/4/5 compliance checks,
#' and organizing allele columns. Automatically detects whether data is from a
#' family study or regular typing format and includes relevant identifiers in the QC output.
#'
#' @param df A data frame or tibble containing raw HLA typing data.
#' @param isfamilydata Logical or NULL. If NULL (default), automatically detects data type.
#'        If TRUE, processes as family data. If FALSE, processes as regular typing data.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A tibble with cleaned and formatted HLA allele data. A QC report is attached
#'         as an attribute (`attr(df, "qc")`) containing sample identifiers and DRB1–DRB3/4/5
#'         compliance flags.
#'
#' @examples
#' \dontrun{
#' raw_data <- read.csv("hla_typing.csv")
#' formatted <- reformat_typing_data(raw_data)
#' qc_report <- attr(formatted, "qc")
#' }
#'
#' @importFrom janitor clean_names
#' @importFrom dplyr across mutate arrange
#' @importFrom forcats fct_relevel
#'
#' @export
reformat_typing_data <- function(df,
                                 isfamilydata = NULL,
                                 quiet = FALSE) {
   # Step 0: Input validation
   if (!is.data.frame(df)) {
      stop("\t❌ Input must be a data frame or tibble.")
   }
   if (nrow(df) == 0) {
      warning("\t⚠️ Input data frame has zero rows.")
   }

   # Step 1: Clean column names
   df <- clean_column_names(df, quiet)

   # Step 2: Standardize allele format
   df <- standardize_allele_format(df, quiet)

   # Step 3: Auto-detect or use user-specified format
   if (is.null(isfamilydata)) {
      if (!quiet) message("\t🔍 Auto-detecting data format...")
      detect_result <- detect_data_type(df, quiet = quiet)
      isfamilydata <- detect_result$is_family
      id_cols <- if (isfamilydata) {
         c(detect_result$family_col, detect_result$member_col)
      } else {
         detect_result$id_cols
      }
   } else {
      if (!quiet) {
         message(
            "\tℹ️ Data format explicitly provided as ",
            ifelse(isfamilydata, "'FAMILY STUDY'", "'REGULAR TYPING'"),
            "; skipping auto-detection."
         )
      }

      # Manually infer ID columns based on user-specified data type
      id_cols <- if (isfamilydata) {
         family_ids <- c("FAMILY_ID", "Family_ID", "family_id", "FamilyID", "familyid", "Family", "family")
         member_ids <- c("Family_Member", "FamilyMember", "Member", "MEMBER", "Relationship", "RELATIONSHIP")
         intersect(c(family_ids, member_ids), names(df))
      } else {
         sample_ids <- c(
            "SampleID", "Sample_ID", "ID", "Id", "id", "PATIENT_ID", "Patient_ID",
            "Subject_ID", "Donor_ID", "donor_id", "RecipientID", "recipient_id"
         )
         intersect(sample_ids, names(df))
      }

      if (length(id_cols) == 0) id_cols <- NULL
   }

   # Step 4: Fill missing second alleles + DRB1-DRB3/4/5 linkage compliance
   fill_result <- fill_missing_alleles(df, id_cols = id_cols, quiet = quiet)
   df <- fill_result$cleaned
   qc <- fill_result$qc

   # Step 5: Organize columns and sort
   df <- organize_and_sort(df, isfamilydata, quiet)

   # Step 6: Attach QC report
   attr(df, "qc") <- qc

   if (!quiet) {
      message("\t✅ Formatting complete. HLA data successfully processed.")
   }

   df
}
