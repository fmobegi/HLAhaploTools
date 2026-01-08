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
