# ~~~~~~~~~~~~~~~~~~~~~~~~~
# standardize_allele_format
# ~~~~~~~~~~~~~~~~~~~~~~~~~
#' Standardize HLA Allele Format
#'
#' Internal helper function to ensure consistent allele notation across columns.
#'
#' @param df A data frame with cleaned column names
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A tibble with standardized allele notations
#' @keywords internal
standardize_allele_format <- function(df, quiet = FALSE) {
  if (!quiet) {
    message("\t🧬 Standardizing HLA allele notation...")
  }

  # Define column types we're expecting
  id_cols <- c("FAMILY_ID", "Family_Member")
  class1 <- c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2")
  class2 <- c(
    "DRB1_1",
    "DRB1_2",
    "DRB3_1",
    "DRB3_2",
    "DRB4_1",
    "DRB4_2",
    "DRB5_1",
    "DRB5_2",
    "DQA1_1",
    "DQA1_2",
    "DQB1_1",
    "DQB1_2",
    "DPA1_1",
    "DPA1_2",
    "DPB1_1",
    "DPB1_2",
    "DMA_1",
    "DMA_2",
    "DMB_1",
    "DMB_2",
    "DOA_1",
    "DOA_2",
    "DOB_1",
    "DOB_2"
  )
  nonclass <- c(
    "F_1",
    "F_2",
    "G_1",
    "G_2",
    "H_1",
    "H_2",
    "J_1",
    "J_2",
    "K_1",
    "K_2",
    "L_1",
    "L_2",
    "E_1",
    "E_2",
    "HFE_1",
    "HFE_2"
  )
  mic <- c("MICA_1", "MICA_2", "MICB_1", "MICB_2")

  detect_group <- function(label, group) {
    present <- intersect(group, names(df))
    if (!quiet) {
      message(
        sprintf(
          if (length(present)) {
            "\t\t✅ %s genes detected: %s"
          } else {
            "\t\t❌ No %s genes detected."
          },
          label,
          paste(unique(sub(
            "_.*", "", present
          )), collapse = ", ")
        )
      )
    }
    present
  }

  # Find allele columns that exist in the dataset
  class1_present <- detect_group("Class I", class1)
  class2_present <- detect_group("Class II", class2)
  nonclass_present <- detect_group("Non-Classical", nonclass)
  mic_present <- detect_group("MIC", mic)

  allele_cols <- c(
    class1_present,
    class2_present,
    nonclass_present,
    mic_present
  )

  if (length(allele_cols) == 0 && !quiet) {
    message("\t\t⚠️ No HLA gene columns detected! Check your column naming.")
  }

  # Standardize allele format for each column
  for (col in allele_cols) {
    gene <- sub("_.*", "", col)
    df[[col]] <- ifelse(
      is.na(df[[col]]) | grepl("^\\s*$", df[[col]]),
      NA_character_,
      ifelse(
        grepl("^\\*", df[[col]]),
        paste0(gene, df[[col]]),
        ifelse(grepl("\\*", df[[col]]), df[[col]], paste0(gene, "*", df[[col]]))
      )
    )
  }

  df
}
