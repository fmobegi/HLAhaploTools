# ~~~~~~~~~~~~~~~~
# load_typing_data
# ~~~~~~~~~~~~~~~~
#' Load HLA Typing Data from CSV, TSV, Excel or RDS
#'
#' Reads an HLA typing data file from various formats and returns a tibble.
#' Interprets common NA-like values (e.g., "", "NULL", "NA", "Unknown") as missing.
#'
#' @param filepath Path to the data file (.csv, .tsv, .txt, .xlsx, or .rds).
#' @param sheet Excel sheet name or index (optional; defaults to the first sheet).
#' @param quiet Logical. If TRUE, suppresses status messages.
#' @param max_file_mb Maximum file size in MB to load without warning (default = 100).
#' @param ... Additional arguments passed to the appropriate reading function.
#'
#' @return A tibble containing the cleaned HLA typing data.
#' @examples
#' \dontrun{
#' # Load from CSV
#' typing_data <- load_typing_data("path/to/hla_data.csv")
#'
#' # Load from Excel with specific sheet
#' excel_data <- load_typing_data("path/to/hla_data.xlsx", sheet = "Family1")
#'
#' # Load silently without messages
#' quiet_load <- load_typing_data("path/to/data.tsv", quiet = TRUE)
#' }
#'
#' @importFrom readr read_csv read_tsv cols
#' @importFrom readxl read_excel
#' @importFrom tibble as_tibble
#' @importFrom tools file_ext
#'
#' @export
load_typing_data <- function(filepath,
                             sheet = NULL,
                             quiet = FALSE,
                             max_file_mb = 100,
                             ...) {
  # Check file existence
  if (!file.exists(filepath)) {
    stop("\t❌ File does not exist: ", filepath)
  }

  # Check file size
  file_size_mb <- file.info(filepath)$size / 1024^2
  if (file_size_mb > max_file_mb) {
    warning(sprintf(
      "\t⚠️ Large file detected (%.1f MB). Loading may take time.",
      file_size_mb
    ))
  }

  ext <- tolower(tools::file_ext(filepath))
  na_values <- c("", NULL, NA, "NULL", "NA", "Unknown")

  if (!quiet) {
    message(sprintf("\t📦 Detected file type: .%s", ext))
  }

  safe_read <- function(read_fun, ...) {
    tryCatch(
      {
        df <- read_fun(..., na = na_values)

        if (!quiet) {
          message(sprintf(
            "\t✅ Loaded file with %d rows and %d columns.",
            nrow(df),
            ncol(df)
          ))
        }

        return(df)
      },
      error = function(e) {
        stop("\t❌ Failed to load file: ", e$message)
      }
    )
  }

  df <- switch(ext,
    "csv" = safe_read(
      readr::read_csv,
      file = filepath,
      col_types = readr::cols(.default = "c"),
      ...
    ),
    "tsv" = safe_read(
      readr::read_tsv,
      file = filepath,
      col_types = readr::cols(.default = "c"),
      ...
    ),
    "txt" = safe_read(
      readr::read_tsv,
      file = filepath,
      col_types = readr::cols(.default = "c"),
      ...
    ),
    "xlsx" = safe_read(
      readxl::read_excel,
      path = filepath,
      sheet = sheet,
      col_types = "text",
      ...
    ),
    "rds" = safe_read(readRDS, file = filepath, ...),
    stop(
      "\t❌ Unsupported file extension: ",
      ext,
      ". Use .csv, .tsv, .txt, .xlsx, or .rds instead."
    )
  )

  df <- tibble::as_tibble(df)
  df
}


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


# ~~~~~~~~~~~~~~~~~~~~
# fill_missing_alleles
# ~~~~~~~~~~~~~~~~~~~~
#' Fill Missing Alleles with DRB1 vs DRB3/4/5 Compliance Checks
#'
#' Internal helper function to:
#' - Copy first allele to second allele position when missing,
#' - Reorder DRB alleles to ensure DRBx_1 is populated when present only in _2,
#' - Fill DRB3/4/5 alleles based on DRB1 linkage rules,
#' - Identify DRB1-DRB3/4/5 inconsistencies.
#'
#' @param df A data frame with standardized allele columns (e.g., A_1, DRB1_2, etc.).
#' @param id_cols Optional character vector of column names to identify samples (e.g., "SampleID").
#'                If NULL, common ID column names will be automatically detected.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{cleaned}{The updated tibble with alleles filled and reordered.}
#'   \item{qc}{A tibble with ID columns and a `DRB145_Compliant` logical flag.}
#' }
#' @keywords internal
fill_missing_alleles <- function(df, id_cols = NULL, quiet = FALSE) {
  if (!quiet) message("\t🔄 Starting allele fill with DRB reordering and compliance check...")

  # Auto-detect ID columns if not provided
  if (is.null(id_cols)) {
    id_candidates <- c(
      "FAMILY_ID", "Family_ID", "family_id", "FamilyID", "familyid", "Family", "family",
      "Family_Member", "FamilyMember", "Member", "MEMBER", "Relationship", "RELATIONSHIP",
      "SampleID", "Sample_ID", "ID", "Id", "id", "PATIENT_ID", "Patient_ID", "Subject_ID",
      "Donor_ID", "donor_id", "RecipientID", "recipient_id"
    )
    id_cols <- intersect(id_candidates, names(df))
  }

  # Helper function to get sample identifier string
  get_sample_id <- function(i) {
    if (is.null(id_cols)) {
      return(paste0("Row", i))
    }
    id_vals <- sapply(id_cols, function(col) {
      val <- df[[col]][i]
      if (is.na(val) || val == "") {
        return(NULL)
      }
      as.character(val)
    })
    id_vals <- id_vals[!sapply(id_vals, is.null)]
    if (length(id_vals) == 0) {
      return(paste0("Row", i))
    }
    paste(id_vals, collapse = "_")
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Step 1: Reorder DRB alleles
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  reorder_drb_alleles <- function(df_local, quiet_local) {
    drb_genes <- c("DRB1", "DRB3", "DRB4", "DRB5")
    swapped_total <- 0

    for (gene in drb_genes) {
      a1 <- paste0(gene, "_1")
      a2 <- paste0(gene, "_2")
      if (all(c(a1, a2) %in% names(df_local))) {
        for (i in seq_len(nrow(df_local))) {
          if ((is.na(df_local[[a1]][i]) || df_local[[a1]][i] == "") &&
            !(is.na(df_local[[a2]][i]) || df_local[[a2]][i] == "")) {
            df_local[[a1]][i] <- df_local[[a2]][i]
            df_local[[a2]][i] <- NA
            swapped_total <- swapped_total + 1
            if (!quiet_local) {
              message(
                sprintf(
                  "\t🔄 Reordered %s alleles for sample %s (swapped _1 and _2)",
                  gene,
                  get_sample_id(i)
                )
              )
            }
          }
        }
      }
    }

    if (!quiet_local && swapped_total == 0) {
      message("\tℹ️ No DRB alleles needed reordering.")
    }
    df_local
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Step 2: Fill missing non-DRB alleles
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  fill_other_genes <- function(df_local, quiet_local) {
    allele_cols <- grep("_[12]$", names(df_local), value = TRUE)
    gene_names <- unique(sub("_[12]$", "", allele_cols))
    non_drb_genes <- setdiff(gene_names, c("DRB1", "DRB3", "DRB4", "DRB5"))

    filled_count <- 0
    for (gene in non_drb_genes) {
      a1 <- paste0(gene, "_1")
      a2 <- paste0(gene, "_2")
      if (all(c(a1, a2) %in% names(df_local))) {
        for (i in seq_len(nrow(df_local))) {
          if ((is.na(df_local[[a2]][i]) || df_local[[a2]][i] == "") &&
            !is.na(df_local[[a1]][i]) && df_local[[a1]][i] != "") {
            df_local[[a2]][i] <- df_local[[a1]][i]
            filled_count <- filled_count + 1
            if (!quiet_local) {
              message(
                sprintf(
                  "\t   Filled missing %s_2 allele for sample %s by copying %s_1",
                  gene,
                  get_sample_id(i),
                  gene
                )
              )
            }
          }
        }
      }
    }

    if (!quiet_local && filled_count == 0) {
      message("\tℹ️ No non-DRB alleles needed filling.")
    }
    df_local
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Step 3: DRB1 vs DRB3/4/5 linkage compliance
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  check_and_fill_drb_linkage <- function(df_local, quiet_local) {
    linkage_map <- list(
      "DRB3" = c("03", "11", "12", "13"),
      "DRB4" = c("07", "08", "09"),
      "DRB5" = c("15", "16")
    )

    extract_group <- function(allele) {
      if (is.na(allele) || allele == "") {
        return(NA)
      }
      m <- regmatches(allele, regexpr("\\*(\\d{2})", allele))
      if (length(m) == 0) {
        return(NA)
      }
      sub("\\*", "", m)
    }

    compliance_vec <- logical(nrow(df_local))
    n_noncompliant <- 0

    for (i in seq_len(nrow(df_local))) {
      drb1_1 <- df_local[["DRB1_1"]][i]
      drb1_2 <- df_local[["DRB1_2"]][i]
      drb1_groups <- na.omit(unique(c(extract_group(drb1_1), extract_group(drb1_2))))
      compliant <- TRUE

      for (drb_gene in names(linkage_map)) {
        a1 <- paste0(drb_gene, "_1")
        a2 <- paste0(drb_gene, "_2")
        if (!all(c(a1, a2) %in% names(df_local))) next

        drb1_requires <- any(drb1_groups %in% linkage_map[[drb_gene]])
        alleles_present <- !all(is.na(df_local[i, c(a1, a2)]) | df_local[i, c(a1, a2)] == "")

        if (drb1_requires && !alleles_present) {
          fill_value <- drb1_1
          if (!is.na(fill_value) && fill_value != "") {
            df_local[i, c(a1, a2)] <- fill_value
            if (!quiet_local) {
              message(
                sprintf(
                  "\t   Filled missing %s alleles for sample %s using DRB1_1 allele",
                  drb_gene, get_sample_id(i)
                )
              )
            }
          } else {
            compliant <- FALSE
            if (!quiet_local) {
              message(
                sprintf(
                  "\t⚠️ Sample %s missing required %s alleles and DRB1 allele not available for fill",
                  get_sample_id(i), drb_gene
                )
              )
            }
          }
        } else if (!drb1_requires && alleles_present) {
          compliant <- FALSE
          if (!quiet_local) {
            message(
              sprintf(
                "\t⚠️ Sample %s has %s alleles present but DRB1 group(s) [%s] does not require it",
                get_sample_id(i), drb_gene, paste(drb1_groups, collapse = "/")
              )
            )
          }
        }
      }

      compliance_vec[i] <- compliant
      if (!compliant) n_noncompliant <- n_noncompliant + 1
    }

    if (!quiet_local && n_noncompliant > 0) {
      message(sprintf("\t⚠️ Found %d samples with DRB1-DRBx linkage inconsistencies.", n_noncompliant))
    }

    list(df = df_local, compliance = compliance_vec)
  }

  # Execute the steps
  df <- reorder_drb_alleles(df, quiet)
  df <- fill_other_genes(df, quiet)
  drb_result <- check_and_fill_drb_linkage(df, quiet)

  # Prepare output
  cleaned_df <- drb_result$df
  compliance_flag <- drb_result$compliance

  qc_df <- if (!is.null(id_cols)) {
    cleaned_df[, id_cols, drop = FALSE]
  } else {
    data.frame(RowNumber = seq_len(nrow(cleaned_df)))
  }
  qc_df$DRB145_Compliant <- compliance_flag

  list(cleaned = cleaned_df, qc = qc_df)
}


# ~~~~~~~~~~~~~~~~~
# organize_and_sort
# ~~~~~~~~~~~~~~~~~
#' Organize and Sort HLA Data
#'
#' Internal helper function to organize columns and sort data.
#' Orders columns by gene class (Class I, Class II, non-classical, MIC)
#' and sorts records by family structure when applicable.
#'
#' @param df A data frame with cleaned and standardized HLA data
#' @param isfamilydata Logical. If TRUE, sort by family member.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A sorted and organized tibble
#' @keywords internal
#' @importFrom dplyr arrange
organize_and_sort <- function(df,
                              isfamilydata = TRUE,
                              quiet = FALSE) {
  if (!quiet) {
    message("\t📋 Organizing columns and sorting data...")
  }

  # Find ID columns based on data type
  id_cols <- if (isfamilydata) {
    # For family data, use family columns
    c("FAMILY_ID", "Family_Member")
  } else {
    # For non-family data, use any ID column we can find
    possible_id_cols <- c(
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
    # If no standard ID column, fallback to FAMILY_ID if present
    intersect(c(possible_id_cols, "FAMILY_ID"), names(df))
  }

  # Find existing allele columns
  allele_cols <- grep("_[12]$", names(df), value = TRUE)

  # Define standard gene orders within each class
  class1_order <- c("A", "B", "C")
  class2_order <- c(
    "DRB1",
    "DRB3",
    "DRB4",
    "DRB5",
    "DQA1",
    "DQB1",
    "DPA1",
    "DPB1",
    "DOA",
    "DOB",
    "DMA",
    "DMB"
  )
  nonclass_order <- c("E", "F", "G", "H", "HFE", "J", "K", "L")
  mic_order <- c("MICA", "MICB")

  # Extract gene names from columns (without _1 or _2 suffix)
  gene_names <- unique(sub("_[12]$", "", allele_cols))

  # Group genes by class
  class1_genes <- intersect(class1_order, gene_names)
  class2_genes <- intersect(class2_order, gene_names)
  nonclass_genes <- intersect(nonclass_order, gene_names)
  mic_genes <- intersect(mic_order, gene_names)
  other_genes <- setdiff(
    gene_names,
    c(class1_order, class2_order, nonclass_order, mic_order)
  )

  # Sort genes within each class according to standard order
  class1_genes <- class1_genes[match(class1_genes, class1_order)]
  class2_genes <- class2_genes[match(class2_genes, class2_order)]
  nonclass_genes <- nonclass_genes[match(nonclass_genes, nonclass_order)]
  mic_genes <- mic_genes[match(mic_genes, mic_order)]
  other_genes <- sort(other_genes) # Alphabetically sort any other genes

  # Combine all genes in order of importance
  ordered_genes <- c(
    class1_genes,
    class2_genes,
    nonclass_genes,
    mic_genes,
    other_genes
  )

  # Get columns in desired order (both _1 and _2 for each gene)
  ordered_allele_cols <- c()
  for (gene in ordered_genes) {
    col1 <- paste0(gene, "_1")
    col2 <- paste0(gene, "_2")
    if (col1 %in% names(df)) {
      ordered_allele_cols <- c(ordered_allele_cols, col1)
    }
    if (col2 %in% names(df)) {
      ordered_allele_cols <- c(ordered_allele_cols, col2)
    }
  }

  # Report gene organization if not quiet
  if (!quiet && length(gene_names) > 0) {
    message("\t   Organizing genes by HLA class:")
    if (length(class1_genes) > 0) {
      message(sprintf("\t   ▪ Class I genes: %s", paste(class1_genes, collapse = ", ")))
    }
    if (length(class2_genes) > 0) {
      message(sprintf("\t   ▪ Class II genes: %s", paste(class2_genes, collapse = ", ")))
    }
    if (length(nonclass_genes) > 0) {
      message(sprintf(
        "\t   ▪ Non-classical genes: %s",
        paste(nonclass_genes, collapse = ", ")
      ))
    }
    if (length(mic_genes) > 0) {
      message(sprintf("\t   ▪ MIC genes: %s", paste(mic_genes, collapse = ", ")))
    }
    if (length(other_genes) > 0) {
      message(sprintf("\t   ▪ Other genes: %s", paste(other_genes, collapse = ", ")))
    }
  }

  # Order columns: ID columns, then ordered allele columns, then everything else
  ordering_cols <- c(intersect(id_cols, names(df)), ordered_allele_cols)
  remaining <- setdiff(names(df), ordering_cols)

  # Apply column ordering
  df <- df[, c(ordering_cols, remaining), drop = FALSE]

  # Remove empty rows and columns
  df <- df[rowSums(is.na(df)) < ncol(df), , drop = FALSE]
  df <- df[, colSums(is.na(df)) < nrow(df), drop = FALSE]

  # Sort by family structure if applicable
  if (isfamilydata &&
    all(c("FAMILY_ID", "Family_Member") %in% names(df))) {
    # Standardize family member codes for sorting
    df$Family_Member <- sub("^[fF]$|^[fF]ather$|^[dD]ad$", "F", df$Family_Member)
    df$Family_Member <- sub("^[mM]$|^[mM]other$|^[mM]om$", "M", df$Family_Member)
    df$Family_Member <- sub("^[cC]hild[ -]?([0-9]+)$", "C\\1", df$Family_Member)
    df$Family_Member <- sub("^[cC]hild$", "C1", df$Family_Member)

    # Standard ordering: Father, Mother, Child1, Child2, ...
    child_codes <- grep("^C[0-9]+$", unique(df$Family_Member), value = TRUE)
    order_levels <- c("F", "M", sort(child_codes))

    df$Family_Member <- factor(df$Family_Member, levels = order_levels)
    df <- dplyr::arrange(df, FAMILY_ID, Family_Member)

    if (!quiet) {
      message("\t👪 Data sorted by family structure.")
    }
  }

  df
}


# ~~~~~~~~~~~~~~~~~~~~
# decode_classical_mac
# ~~~~~~~~~~~~~~~~~~~~
#' Decode MAC-Encoded Alleles in Classical HLA Columns
#'
#' Uses immunotation to decode multi-allele codes (MAC) in classical HLA loci.
#'
#' @param df A data frame with HLA allele columns (e.g. A_1, B_2, etc).
#' @param loci Optional set of classical genes to decode.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return Data frame with decoded MAC values.
#' @examples
#' \dontrun{
#' # Decode all classical loci
#' typed_data <- load_typing_data("hla_typing.csv") %>%
#'   reformat_typing_data() %>%
#'   decode_classical_mac()
#'
#' # Only decode Class I loci
#' class1_decoded <- decode_classical_mac(typed_data, loci = c("A", "B", "C"))
#' }
#' @export
#' @importFrom immunotation decode_MAC
decode_classical_mac <- function(df,
                                 loci = c(
                                   "A",
                                   "B",
                                   "C",
                                   "DRB1",
                                   "DRB3",
                                   "DRB4",
                                   "DRB5",
                                   "DQA1",
                                   "DQB1",
                                   "DPA1",
                                   "DPB1"
                                 ),
                                 quiet = FALSE) {
  # Check for immunotation package
  if (!requireNamespace("immunotation", quietly = TRUE)) {
    stop(
      "\t❌ The 'immunotation' package is required. Please install it with BiocManager::install('immunotation')."
    )
  }

  if (!quiet) {
    message("\t🔍 Decoding MAC strings in HLA alleles...")
  }

  # Find columns to process
  allele_cols <- intersect(unlist(lapply(loci, function(l) {
    paste0(l, c(
      "_1", "_2"
    ))
  })), names(df))

  if (length(allele_cols) == 0) {
    warning("\t⚠️ No matching HLA columns found for MAC decoding.")
    return(df)
  }

  # MAC pattern
  mac_pattern <- "^[A-Z0-9]+\\*[0-9]{2}:[A-Z]+$"
  mac_count <- 0

  # Process each column
  for (col in allele_cols) {
    df[[col]] <- sapply(df[[col]], function(allele) {
      if (is.na(allele) || allele == "") {
        return(allele)
      }
      if (grepl(mac_pattern, allele)) {
        mac_count <<- mac_count + 1
        tryCatch(
          immunotation::decode_MAC(allele),
          error = function(e) {
            warning(sprintf("\t⚠️ Failed to decode MAC: %s", allele))
            allele
          }
        )
      } else {
        allele
      }
    })
  }

  if (!quiet) {
    if (mac_count > 0) {
      message(sprintf("\t✅ Decoded %d MAC strings.", mac_count))
    } else {
      message("\tℹ️ No MAC strings found to decode.")
    }
  }

  df
}


# ~~~~~~~~~~~~~~~~
# trim_hla_results
# ~~~~~~~~~~~~~~~~
#' Trim Resolution of HLA Alleles Across a Data Frame
#'
#' Applies `HLAtools::multiAlleleTrim()` to all columns in a tibble or data frame,
#' reducing allele resolution (e.g., from 4-field to 2-field notation).
#'
#' Handles MAC-decoded alleles containing "/" separators by trimming each component.
#' NA values and empty strings are preserved, and trimming failures are gracefully handled with fallback warnings.
#'
#' @param df A tibble or data frame with HLA alleles as character columns.
#' @param resolution Integer. Desired resolution level (e.g., 2 or 4). Default is 3.
#' @param append Logical. If TRUE (default), ambiguity suffixes (e.g., "G", "N") are retained after trimming.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A tibble with trimmed allele values.
#'
#' @examples
#' \dontrun{
#' # Load and trim alleles to 2-field resolution
#' typed_data <- load_typing_data("hla_typing.csv") %>%
#'   reformat_typing_data() %>%
#'   trim_hla_results(resolution = 2)
#' }
#' @importFrom HLAtools multiAlleleTrim
#' @importFrom tibble as_tibble
#'
#' @export
trim_hla_results <- function(df,
                             resolution = 3,
                             append = TRUE,
                             quiet = FALSE) {
  # Check for HLAtools package
  if (!requireNamespace("HLAtools", quietly = TRUE)) {
    stop(
      "\t❌ The 'HLAtools' package is required. Please install it from GitHub: remotes::install_github('DKMS/HLAtools')"
    )
  }

  if (!quiet) {
    message(sprintf("\t✂️ Trimming alleles to %d-field resolution...", resolution))
  }

  # Count of trimmed alleles for reporting
  trim_count <- 0
  fail_count <- 0
  mac_count <- 0

  # Apply trimming to each column
  df[] <- lapply(df, function(col) {
    if (!is.character(col)) {
      return(col)
    }

    sapply(col, function(allele) {
      if (is.na(allele) || allele == "") {
        return(allele)
      }

      # Check if this is a MAC-decoded allele containing "/"
      if (grepl("/", allele)) {
        mac_count <<- mac_count + 1

        # Split by "/" and trim each component
        components <- unlist(strsplit(allele, "/"))
        trimmed_components <- sapply(components, function(component) {
          tryCatch(
            {
              trimmed <- HLAtools::multiAlleleTrim(component,
                resolution = resolution,
                append = append
              )
              trim_count <<- trim_count + 1
              return(trimmed)
            },
            error = function(e) {
              fail_count <<- fail_count + 1
              if (!quiet) {
                message(
                  sprintf(
                    "\t⚠️ Trimming failed for '%s' in MAC-decoded string — keeping original.",
                    component
                  )
                )
              }
              return(component)
            }
          )
        })

        # Recombine with "/"
        return(paste(trimmed_components, collapse = "/"))
      }

      # Regular allele (no MAC decoding)
      tryCatch(
        {
          trimmed <- HLAtools::multiAlleleTrim(allele, resolution = resolution, append = append)
          trim_count <<- trim_count + 1
          return(trimmed)
        },
        error = function(e) {
          fail_count <<- fail_count + 1
          if (!quiet) {
            message(sprintf(
              "\t⚠️ Trimming failed for '%s' — keeping original.",
              allele
            ))
          }
          return(allele)
        }
      )
    })
  })

  if (!quiet) {
    message(
      sprintf(
        "\t✅ Trimmed %d alleles to %d-field resolution.",
        trim_count - fail_count,
        resolution
      )
    )
    if (mac_count > 0) {
      message(
        sprintf(
          "\t🧩 Processed %d MAC-decoded alleles with multiple options.",
          mac_count
        )
      )
    }
    if (fail_count > 0) {
      message(sprintf("\t⚠️ %d alleles could not be trimmed properly.", fail_count))
    }
  }
  df <- tibble::as_tibble(df)
  df
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# validate_family_data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Validate Family Data Structure
#'
#' Checks whether the input dataset reflects a valid family study format and
#' provides structural insights. This includes detecting father, mother, and child roles
#' based on the `Family_Member` column and reporting completeness of family trios.
#'
#' @param df A data frame or tibble containing HLA typing data.
#' @param stop_if_invalid Logical. If TRUE, stops execution if data is not valid family format.
#' @param verbose Logical. If TRUE, prints detailed validation output.
#'
#' @return Logical `TRUE` if valid family structure is detected, `FALSE` otherwise.
#' @keywords internal
validate_family_data <- function(df,
                                 stop_if_invalid = TRUE,
                                 verbose = TRUE) {
  # Detect data type using quiet mode
  detect <- detect_data_type(df, quiet = !verbose)
  is_family <- detect$is_family

  if (!is_family) {
    if (stop_if_invalid) {
      stop(
        "\n\t❌ Error: Data does not appear to be in family study format.\n",
        "\t   Use other functions for registry or population data."
      )
    }
    return(FALSE)
  }

  if (verbose) {
    message("\n\t🔍 Validating family structure...")

    family_col <- "FAMILY_ID"
    member_col <- "Family_Member"

    # Summary counts
    n_families <- length(unique(df[[family_col]]))
    n_fathers <- sum(df[[member_col]] %in% c("F", "Father"), na.rm = TRUE)
    n_mothers <- sum(df[[member_col]] %in% c("M", "Mother"), na.rm = TRUE)
    n_children <- sum(grepl("^C\\d+$|Child", df[[member_col]]), na.rm = TRUE)

    message(sprintf(
      "\t   Found %d families with %d fathers, %d mothers, %d children.",
      n_families, n_fathers, n_mothers, n_children
    ))

    # Trio completeness check
    complete_families <- sum(sapply(unique(df[[family_col]]), function(fam) {
      fam_data <- df[df[[family_col]] == fam, ]
      has_father <- any(fam_data[[member_col]] %in% c("F", "Father"))
      has_mother <- any(fam_data[[member_col]] %in% c("M", "Mother"))
      has_child <- any(grepl("^C\\d+$|Child", fam_data[[member_col]]))
      has_father && has_mother && has_child
    }))

    message(sprintf(
      "\t   %d/%d families have complete trios (father, mother, child).",
      complete_families, n_families
    ))
  }

  TRUE
}


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


# ~~~~~~~~~~~~~~~~~~~~~
# check_deleted_alleles
# ~~~~~~~~~~~~~~~~~~~~~
#' Check HLA dataframe for deleted alleles
#'
#' @description Analyzes a dataframe containing HLA allele designations, checking for
#' deleted alleles according to the IMGT/HLA database. Downloads and parses the official
#' deleted allele list, then matches any occurrences in the user's data.
#'
#' @param df_hla A dataframe containing HLA allele typing data.
#' @param deleted_file_path Path or URL to the Deleted_alleles.txt file
#'        (default: IMGT/HLA database URL).
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A dataframe containing information about deleted alleles found in the input.
#'
#' @examples
#' \dontrun{
#' results <- check_deleted_alleles(df_hla)
#' }
#' @keywords internal
#' @importFrom curl curl
check_deleted_alleles <- function(
    df_hla,
    deleted_file_path = "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Deleted_alleles.txt",
    quiet = FALSE) {
  if (!quiet) {
    message("\t🔍 Checking for deleted HLA alleles...")
    message("\t📥 Downloading current list of deleted alleles from:\n\t ", deleted_file_path)
  }

  deleted_alleles <- list()

  tryCatch(
    {
      con <- if (grepl("^https?://", deleted_file_path)) {
        curl::curl(deleted_file_path)
      } else {
        file(deleted_file_path, open = "r")
      }

      lines <- readLines(con, warn = FALSE)
      close(con)

      header_found <- FALSE

      for (line in lines) {
        if (grepl("^#", line) || nchar(trimws(line)) == 0) next
        if (grepl("^AlleleID,Allele,Description", line)) {
          header_found <- TRUE
          next
        }
        if (header_found) {
          split_line <- strsplit(line, ",")[[1]]
          if (length(split_line) >= 2) {
            allele_name <- trimws(split_line[2])
            description <- if (length(split_line) >= 3) {
              paste(split_line[3:length(split_line)], collapse = ",")
            } else {
              "No description"
            }
            deleted_alleles[[allele_name]] <- description
          }
        }
      }

      if (!quiet) {
        message(paste0("\t📊 Loaded ", length(deleted_alleles), " deleted alleles."))
      }
    },
    error = function(e) {
      message("\t❌ Error reading deleted alleles file: ", e$message)
      data.frame()
    }
  )

  if (length(deleted_alleles) == 0) {
    if (!quiet) {
      message("\t⚠️ No deleted alleles found in the reference file.")
    }
    data.frame()
  }

  results <- data.frame(
    FAMILY_ID = character(),
    Family_Member = character(),
    Column = character(),
    Deleted_Allele = character(),
    Description = character(),
    stringsAsFactors = FALSE
  )

  hla_columns <- setdiff(names(df_hla), c("FAMILY_ID", "Family_Member"))

  for (i in seq_len(nrow(df_hla))) {
    for (col in hla_columns) {
      allele <- df_hla[i, col]
      if (is.na(allele) || allele == "") next
      if (allele %in% names(deleted_alleles)) {
        results <- rbind(results, data.frame(
          FAMILY_ID = df_hla$FAMILY_ID[i],
          Family_Member = df_hla$Family_Member[i],
          Column = col,
          Deleted_Allele = allele,
          Description = deleted_alleles[[allele]],
          stringsAsFactors = FALSE
        ))

        if (!quiet) {
          message(sprintf(
            "\t🔴 Found deleted allele: %s in FAMILY_ID=%s, Family_Member=%s, Column=%s",
            allele, df_hla$FAMILY_ID[i], df_hla$Family_Member[i], col
          ))
        }
      }
    }
  }

  if (!quiet) {
    if (nrow(results) == 0) {
      message("\t✅ No deleted alleles found in the data.")
    } else {
      message(sprintf("\t⚠️ Found %d instances of deleted alleles.", nrow(results)))
    }
  }

  results
}
