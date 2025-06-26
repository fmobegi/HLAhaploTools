# ~~~~~~~~~~~~~~~~~~~~
# load_typing_data()
# ~~~~~~~~~~~~~~~~~~~~
#' Load HLA Typing Data from CSV, TSV, or Excel
#'
#' Reads an HLA typing data file from CSV, TSV, or Excel format and returns a tibble.
#' Interprets common NA-like values (e.g., "", "NULL", "NA", "Unknown") as missing.
#'
#' @param filepath Path to the data file (.csv, .tsv, .txt, or .xlsx).
#' @param sheet Excel sheet name or index (optional; defaults to the first sheet).
#' @param ... Additional arguments passed to the appropriate reading function.
#'
#' @return A tibble containing the cleaned HLA typing data.
#' @importFrom readr read_csv read_tsv cols
#' @importFrom readxl read_excel
#' @importFrom tibble as_tibble
#' @export
load_typing_data <- function(filepath, sheet = NULL, ...) {
  if (!file.exists(filepath)) stop("\t❌ File does not exist: ", filepath)

  ext <- tolower(tools::file_ext(filepath))
  na_values <- c("", NULL, NA, "NULL", "NA", "Unknown")
  message(sprintf("\t📦 Detected file type: .%s", ext))

  safe_read <- function(read_fun, ...) {
    tryCatch(
      {
        df <- read_fun(..., na = na_values)
        message(sprintf("\t✅ Loaded file with %d rows and %d columns.", nrow(df), ncol(df)))
        return(df)
      },
      error = function(e) {
        stop("\t❌ Failed to load file: ", e$message)
      }
    )
  }

  df <- switch(ext,
    "csv"  = safe_read(readr::read_csv, file = filepath, col_types = readr::cols(.default = "c"), ...),
    "tsv"  = safe_read(readr::read_tsv, file = filepath, col_types = readr::cols(.default = "c"), ...),
    "txt"  = safe_read(readr::read_tsv, file = filepath, col_types = readr::cols(.default = "c"), ...),
    "xlsx" = safe_read(readxl::read_excel, path = filepath, sheet = sheet, col_types = "text", ...),
    stop("\t❌ Unsupported file extension: ", ext, ". Use .csv, .tsv, .txt, or .xlsx instead.")
  )

  tibble::as_tibble(df)
}

# ~~~~~~~~~~~~~~~~~~~~
# reformat_typing_data()
# ~~~~~~~~~~~~~~~~~~~~
#' Reformat and Verbosely Validate Raw HLA Typing Data
#'
#' Cleans and formats HLA typing data by standardizing allele notations,
#' filling missing second allele calls, and reordering records.
#'
#' @param df A data frame or tibble containing raw HLA typing data.
#' @param isFamilyData Logical. Default TRUE. If FALSE, Family_Member sorting is skipped.
#'
#' @return A tibble with cleaned and reported HLA allele fields.
#' @export
#' @importFrom janitor clean_names
#' @importFrom dplyr across mutate arrange
#' @importFrom forcats fct_relevel
reformat_typing_data <- function(df, isFamilyData = TRUE) {
  if (!is.data.frame(df)) stop("\t❌ Input must be a data frame or tibble.")

  df <- janitor::clean_names(df, case = "none")
  df <- dplyr::mutate(df, dplyr::across(where(is.character), trimws))

  id_cols <- if (isFamilyData) c("FAMILY_ID", "Family_Member") else "FAMILY_ID"

  class1 <- c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2")
  class2 <- c(
    "DRB1_1", "DRB1_2", "DRB3_1", "DRB3_2", "DRB4_1", "DRB4_2",
    "DRB5_1", "DRB5_2", "DQA1_1", "DQA1_2", "DQB1_1", "DQB1_2",
    "DPA1_1", "DPA1_2", "DPB1_1", "DPB1_2"
  )
  nonclass <- c("F_1", "F_2", "G_1", "G_2", "H_1", "H_2", "J_1", "J_2", "E_1", "E_2")
  mic <- c("MICA_1", "MICA_2", "MICB_1", "MICB_2")

  detect_group <- function(label, group) {
    present <- intersect(group, names(df))
    message(sprintf(
      if (length(present)) "\t✅ %s genes detected: %s" else "\t❌ No %s genes detected.",
      label, paste(unique(sub("_.*", "", present)), collapse = ", ")
    ))
    present
  }

  class1_present <- detect_group("Class I", class1)
  class2_present <- detect_group("Class II", class2)
  nonclass_present <- detect_group("Non-Classical", nonclass)
  mic_present <- detect_group("MIC", mic)

  allele_cols <- c(class1_present, class2_present, nonclass_present, mic_present)
  ordering_cols <- c(intersect(id_cols, names(df)), allele_cols)
  remaining <- setdiff(names(df), ordering_cols)

  for (col in allele_cols) {
    gene <- sub("_.*", "", col)
    df[[col]] <- ifelse(is.na(df[[col]]) | grepl("^\\s*$", df[[col]]), NA_character_,
      ifelse(grepl("^\\*", df[[col]]), paste0(gene, df[[col]]),
        ifelse(grepl("\\*", df[[col]]), df[[col]], paste0(gene, "*", df[[col]]))
      )
    )
  }

  for (gene in unique(sub("_.*", "", allele_cols))) {
    a1 <- paste0(gene, "_1")
    a2 <- paste0(gene, "_2")
    if (a1 %in% names(df) && a2 %in% names(df)) {
      df[[a2]] <- ifelse((is.na(df[[a2]]) | df[[a2]] == "") & !is.na(df[[a1]]) & df[[a1]] != "", df[[a1]], df[[a2]])
    }
  }

  df <- df[, c(ordering_cols, remaining), drop = FALSE]
  df <- df[rowSums(is.na(df)) < ncol(df), , drop = FALSE]
  df <- df[, colSums(is.na(df)) < nrow(df), drop = FALSE]

  if (isFamilyData && all(c("FAMILY_ID", "Family_Member") %in% names(df))) {
    df <- dplyr::mutate(df, Family_Member = forcats::fct_relevel(Family_Member, c("F", "M", "C1", "C2")))
    df <- dplyr::arrange(df, FAMILY_ID, Family_Member)
  }

  message("\t✅ Formatting complete. HLA data successfully processed.")
  return(df)
}

# ~~~~~~~~~~~~~~~~~~~~
# decode_classical_mac()
# ~~~~~~~~~~~~~~~~~~~~
#' Decode MAC-Encoded Alleles in Classical HLA Columns
#'
#' Uses immunotation to decode multi-allele codes (MAC) in classical HLA loci.
#'
#' @param df A data frame with HLA allele columns (e.g. A_1, B_2, etc).
#' @param loci Optional set of classical genes to decode.
#'
#' @return Data frame with decoded MAC values.
#' @export
#' @importFrom immunotation decode_MAC
decode_classical_mac <- function(df,
                                 loci = c("A", "B", "C", "DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1")) {
  if (!requireNamespace("immunotation", quietly = TRUE)) {
    stop("The 'immunotation' package is required. Please install it with BiocManager::install('immunotation').")
  }

  allele_cols <- intersect(unlist(lapply(loci, function(l) paste0(l, c("_1", "_2")))), names(df))
  mac_pattern <- "^[A-Z0-9]+\\*[0-9]{2}:[A-Z]+$"

  for (col in allele_cols) {
    df[[col]] <- sapply(df[[col]], function(allele) {
      if (is.na(allele) || allele == "") {
        return(allele)
      }
      if (grepl(mac_pattern, allele)) {
        tryCatch(immunotation::decode_MAC(allele),
          error = function(e) {
            warning(sprintf("⚠️ Failed to decode MAC: %s", allele))
            allele
          }
        )
      } else {
        allele
      }
    })
  }
  return(df)
}

# ~~~~~~~~~~~~~~~~~~~~
# trim_hla_results()
# ~~~~~~~~~~~~~~~~~~~~
#' Trim Resolution of HLA Alleles Across a Data Frame
#'
#' Applies `HLAtools::multiAlleleTrim()` to all columns in a tibble or data frame,
#' reducing allele resolution (e.g., from 4-field to 2-field notation).
#'
#' NA values and empty strings are preserved, and trimming failures are gracefully handled with fallback warnings.
#'
#' @param df A tibble or data frame with HLA alleles as character columns.
#' @param resolution Integer. Desired resolution level (e.g., 2 or 4). Default is 3.
#' @param append Logical. If TRUE (default), ambiguity suffixes (e.g., "G", "N") are retained after trimming.
#'
#' @return A tibble with trimmed allele values.
#' @export
#' @importFrom HLAtools multiAlleleTrim
#' @importFrom tibble as_tibble
trim_hla_results <- function(df, resolution = 3, append = TRUE) {
  df[] <- lapply(df, function(col) {
    sapply(col, function(allele) {
      if (is.na(allele) || allele == "") {
        return(allele)
      }
      tryCatch(
        HLAtools::multiAlleleTrim(allele, resolution = resolution, append = append),
        error = function(e) {
          message(sprintf("\t⚠️ Trimming failed for '%s' — keeping original.", allele))
          allele
        }
      )
    })
  })

  message("\t✅ Allele resolution trimming complete.")
  return(tibble::as_tibble(df))
}
