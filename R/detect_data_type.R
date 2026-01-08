# ~~~~~~~~~~~~~~~~
# detect_data_type
# ~~~~~~~~~~~~~~~~
#' Detect HLA Data Type
#'
#' Determines whether data appears to be from a family study or regular typing.
#'
#' @param df A data frame with HLA typing data
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return Logical TRUE if data appears to be family data, FALSE otherwise
#' @keywords internal
detect_data_type <- function(df, quiet = FALSE) {
  col_names <- colnames(df)
  col_names_lower <- tolower(col_names)

  # Define column name patterns (lowercase only for matching)
  family_id_keys <- c("family_id", "familyid", "family")
  member_keys <- c("family_member", "member", "relationship")
  sample_id_keys <- c("sampleid", "sample_id", "id", "patient_id", "subject_id", "donor_id", "recipientid")

  # Identify matching columns
  family_matches <- which(col_names_lower %in% family_id_keys)
  member_matches <- which(col_names_lower %in% member_keys)
  sample_id_matches <- which(col_names_lower %in% sample_id_keys)

  family_col <- if (length(family_matches) > 0) col_names[family_matches[1]] else NULL
  member_col <- if (length(member_matches) > 0) col_names[member_matches[1]] else NULL
  sample_id_cols <- if (length(sample_id_matches) > 0) col_names[sample_id_matches] else character(0)

  # Strictly check if data qualifies as family format
  is_family <- FALSE
  if (!is.null(family_col) && !is.null(member_col) &&
    family_col %in% colnames(df) && member_col %in% colnames(df)) {
    fam_values <- df[[family_col]]
    member_vals <- unique(na.omit(as.character(df[[member_col]])))
    fam_repeats <- any(table(fam_values) > 1)

    valid_member_labels <- sum(member_vals %in% c("F", "M", "Father", "Mother", "Child")) +
      sum(grepl("^C\\d+$", member_vals))

    is_family <- fam_repeats && valid_member_labels >= 2
  }

  # Fallback logic for non-family format
  id_cols <- if (is_family) {
    na.omit(c(family_col, member_col))
  } else if (length(sample_id_cols) > 0) {
    sample_id_cols
  } else if (!is.null(family_col) && family_col %in% colnames(df)) {
    family_col
  } else {
    fallback_col <- "Sample_ID"
    if (!(fallback_col %in% colnames(df))) {
      df[[fallback_col]] <- seq_len(nrow(df)) # inject row-based IDs
      if (!quiet) message(sprintf("\t⚠️ No ID columns found; added fallback ID column '%s'", fallback_col))
    }
    fallback_col
  }

  # Messages
  if (!quiet) {
    if (is_family) {
      message("\t👪 Data appears to be FAMILY STUDY format")
      message(sprintf("\t   Found family ID column: %s", family_col))
      message(sprintf("\t   Found family member column: %s", member_col))
    } else {
      message("\t👤 Data appears to be REGULAR TYPING format (e.g., registry or population data)")
      if (length(sample_id_cols) > 0) {
        message(sprintf("\t   Found sample ID columns: %s", paste(sample_id_cols, collapse = ", ")))
      } else if (!is.null(family_col)) {
        message(sprintf("\t   Using %s as sample identifier", family_col))
      } else {
        message("\t   No standard identifier column found - will use row numbers")
      }
    }
  }

  # Return metadata list
  list(
    is_family = is_family,
    family_col = if (is_family) family_col else NULL,
    member_col = if (is_family) member_col else NULL,
    id_cols = id_cols
  )
}
