# ~~~~~~~~~~~~~~~~~~
# clean_column_names
# ~~~~~~~~~~~~~~~~~~
#' Clean HLA Data Column Names
#'
#' Internal helper function to standardize column names.
#'
#' @param df A data frame or tibble
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A tibble with cleaned column names
#' @keywords internal
#'
#' @importFrom janitor clean_names
#' @importFrom dplyr mutate across
clean_column_names <- function(df, quiet = FALSE) {
  if (!quiet) {
    message("\t🧹 Cleaning column names...")
  }

  # Make backup of original column names for reporting
  orig_names <- colnames(df)

  # Trim whitespace from character columns
  df <- dplyr::mutate(df, dplyr::across(where(is.character), trimws))

  # Define known HLA gene patterns
  hla_gene_patterns <- c(
    # Class I
    "^[A-C]$",
    "^[A-C][_.-]?\\d$",
    "^[A-C][_.-]?ALLELE[_.-]?\\d$",
    # Class II
    "^DR[A-Z]?\\d*$",
    "^DR[A-Z]?\\d*[_.-]?\\d$",
    "^DR[A-Z]?\\d*[_.-]?ALLELE[_.-]?\\d$",
    "^D[QP][A-Z]\\d*$",
    "^D[QP][A-Z]\\d*[_.-]?\\d$",
    "^D[QP][A-Z]\\d*[_.-]?ALLELE[_.-]?\\d$",
    # Non-classical
    "^[E-L]$",
    "^[E-L][_.-]?\\d$",
    "^[E-L][_.-]?ALLELE[_.-]?\\d$",
    "^HFE$",
    "^HFE[_.-]?\\d$",
    "^HFE[_.-]?ALLELE[_.-]?\\d$",
    # MIC
    "^MIC[AB]$",
    "^MIC[AB][_.-]?\\d$",
    "^MIC[AB][_.-]?ALLELE[_.-]?\\d$"
  )

  # First pass with janitor but preserving case
  df <- janitor::clean_names(df, case = "none")

  # Track changes
  changes <- vector("character", 0)

  # For each column, check if it matches an HLA gene pattern
  for (i in seq_along(colnames(df))) {
    col_name <- colnames(df)[i]

    # Skip columns already in proper format (with _1 or _2 suffix)
    if (grepl("^\\w+_[12]$", col_name)) {
      next
    }

    # Check if column matches HLA gene pattern
    is_hla_gene <- any(sapply(hla_gene_patterns, function(pattern) {
      grepl(pattern, col_name, ignore.case = TRUE)
    }))

    if (is_hla_gene) {
      # Extract base gene name and allele number
      if (grepl("ALLELE[_.-]?(\\d)$", col_name, ignore.case = TRUE)) {
        # Format: GENE_ALLELE_1 or GENE.ALLELE.1
        new_name <- sub("(\\w+)[_.-]?ALLELE[_.-]?(\\d)$",
          "\\1_\\2",
          col_name,
          ignore.case = TRUE
        )
      } else if (grepl("(\\w+)[_.-]?(\\d)$", col_name)) {
        # Format: GENE_1 or GENE.1 or GENE1
        new_name <- sub("(\\w+)[_.-]?(\\d)$", "\\1_\\2", col_name)
      } else if (grepl("^([A-Z]+\\d*)$", col_name, ignore.case = TRUE)) {
        # Format: GENE with no number
        # In this case, we're not sure if this is a gene or not
        # Don't change it
        next
      } else {
        # Can't determine format
        next
      }

      # Record change if any
      if (new_name != col_name) {
        changes <- c(changes, sprintf("'%s' → '%s'", col_name, new_name))
        colnames(df)[i] <- new_name
      }
    }
  }

  # Standardize ID columns

  # (Family ID)
  family_id_variants <- c(
    "FAMILY_ID",
    "Family_ID",
    "family_id",
    "FamilyID",
    "familyid",
    "Family",
    "family"
  )
  for (variant in family_id_variants) {
    idx <- which(tolower(colnames(df)) == tolower(variant))
    if (length(idx) > 0) {
      old_name <- colnames(df)[idx[1]]
      colnames(df)[idx[1]] <- "FAMILY_ID"
      if (old_name != "FAMILY_ID") {
        changes <- c(changes, sprintf("'%s' → 'FAMILY_ID'", old_name))
      }
      break
    }
  }

  # (Family Member)
  member_variants <- c(
    "Family_Member",
    "FamilyMember",
    "Member",
    "MEMBER",
    "Relationship",
    "RELATIONSHIP"
  )
  for (variant in member_variants) {
    idx <- which(tolower(colnames(df)) == tolower(variant))
    if (length(idx) > 0) {
      old_name <- colnames(df)[idx[1]]
      colnames(df)[idx[1]] <- "Family_Member"
      if (old_name != "Family_Member") {
        changes <- c(changes, sprintf("'%s' → 'Family_Member'", old_name))
      }
      break
    }
  }

  # (Sample ID - for non-family data)
  sample_id_variants <- c(
    "SampleID",
    "Sample_ID",
    "ID",
    "Id",
    "id",
    "PATIENT_ID",
    "Patient_ID",
    "Subject_ID",
    "Donor_ID",
    "donor_id",
    "RecipientID",
    "recipient_id"
  )

  for (variant in sample_id_variants) {
    idx <- which(tolower(colnames(df)) == tolower(variant))
    if (length(idx) > 0) {
      old_name <- colnames(df)[idx[1]]
      colnames(df)[idx[1]] <- "Sample_ID"
      if (old_name != "Sample_ID") {
        changes <- c(changes, sprintf("'%s' → 'Sample_ID'", old_name))
      }
      break
    }
  }

  # Report changes
  if (!quiet && length(changes) > 0) {
    message(sprintf("\t   Standardized %d column names", length(changes)))
    if (length(changes) <= 5) {
      for (change in changes) {
        message(sprintf("\t   %s", change))
      }
    } else {
      message(sprintf("\t   (First 5 changes shown of %d total)", length(changes)))
      for (change in changes[1:5]) {
        message(sprintf("\t   %s", change))
      }
    }
  }

  df
}
