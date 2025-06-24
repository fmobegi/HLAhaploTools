#' Load HLA Typing Data from CSV or TSV Using readr
#'
#' Reads an HLA typing data file in CSV or TSV format and returns a tibble.
#'
#' @param filepath Path to the data file (.csv or .tsv).
#' @param ... Additional arguments passed to the readr reading functions.
#'
#' @return A tibble containing the HLA typing data.
#' @importFrom readr read_csv read_tsv
#' @export
load_typing_data <- function(filepath, ...) {
  if (!file.exists(filepath)) {
    stop("File does not exist: ", filepath)
  }

  ext <- tolower(tools::file_ext(filepath))

  # Helper to load file safely
  safe_read <- function(read_fun, path, ...) {
    tryCatch({
      df <- read_fun(path, ...)
      cat("Successfully loaded data with",
          nrow(df),
          "rows and",
          ncol(df),
          "columns\n")
      df
    }, error = function(e) {
      stop("Failed to load file '", path, "' due to error: ", e$message)
    })
  }

  if (ext == "csv") {
    df <- safe_read(readr::read_csv, filepath, ...)
  } else if (ext %in% c("tsv", "txt")) {
    df <- safe_read(readr::read_tsv, filepath, ...)
  } else {
    stop("Unsupported file extension: ",
         ext,
         ". Please use .csv or .tsv files.")
  }

  return(df)
}
