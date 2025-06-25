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
  if (!file.exists(filepath)) {
    stop("File does not exist: ", filepath)
  }

  ext <- tolower(tools::file_ext(filepath))
  na_values <- c("", NULL, NA, "NULL", "NA", "Unknown")

  safe_read <- function(read_fun, ...) {
    tryCatch(
      {
        df <- read_fun(..., na = na_values)
        cat("Successfully loaded data with", nrow(df), "rows and", ncol(df), "columns\n")
        df
      },
      error = function(e) {
        stop("Failed to load file due to error: ", e$message)
      }
    )
  }

  if (ext == "csv") {
    df <- safe_read(readr::read_csv, file = filepath, col_types = readr::cols(.default = "c"), ...)
  } else if (ext %in% c("tsv", "txt")) {
    df <- safe_read(readr::read_tsv, file = filepath, col_types = readr::cols(.default = "c"), ...)
  } else if (ext == "xlsx") {
    df <- safe_read(readxl::read_excel, path = filepath, sheet = sheet, col_types = "text", ...)
  } else {
    stop("Unsupported file extension: ", ext, ". Use .csv, .tsv, .txt, or .xlsx instead.")
  }

  return(tibble::as_tibble(df))
}
