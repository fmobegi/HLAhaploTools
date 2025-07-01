
# Test preprocessing functions with real data files

# Find all data files
data_dir <- system.file("data", package = "HLAhaploTools")
if (data_dir == "") {
  # If not installed, try local path
  data_dir <- "data"
}

# Test data loading from different formats
test_that("load_typing_data works with different file formats", {
  # Try to load each format
  formats <- c("csv", "tsv", "txt", "xlsx")

  for (fmt in formats) {
    file_path <- file.path(data_dir, paste0("family_typing_data.", fmt))

    # Skip if file doesn't exist
    if (!file.exists(file_path)) {
      message("Skipping ", fmt, " test (file not found)")
      next
    }

    # Test loading
    result <- load_typing_data(file_path, quiet = TRUE)
    expect_true(is.data.frame(result))
    expect_gt(nrow(result), 0)
    expect_gt(ncol(result), 0)

    # Check if expected columns exist (adjust based on your data)
    expect_true(any(grepl("_[12]$", names(result))))
  }
})

test_that("standardize_colnames works with real data", {
  # Load a data file (any format)
  formats <- c("csv", "tsv", "txt", "xlsx")
  loaded <- FALSE

  for (fmt in formats) {
    file_path <- file.path(data_dir, paste0("family_typing_data.", fmt))
    if (file.exists(file_path)) {
      raw_data <- load_typing_data(file_path, quiet = TRUE)
      loaded <- TRUE
      break
    }
  }

  if (!loaded) {
    skip("No data files found for testing")
  }

  # Test standardization
  std_data <- standardize_colnames(raw_data, quiet = TRUE)

  # Verify column names are standardized
  expect_false(any(grepl("\\.", names(std_data))))
  expect_false(any(grepl("-", names(std_data))))
  expect_false(any(grepl(" ", names(std_data))))

  # Check for expected HLA columns
  hla_cols <- sum(grepl("_[12]$", names(std_data)))
  expect_gt(hla_cols, 0)
})

test_that("extract_loci works with real data", {
  # Load and standardize a data file
  formats <- c("csv", "tsv", "txt", "xlsx")
  loaded <- FALSE

  for (fmt in formats) {
    file_path <- file.path(data_dir, paste0("family_typing_data.", fmt))
    if (file.exists(file_path)) {
      raw_data <- load_typing_data(file_path, quiet = TRUE)
      std_data <- standardize_colnames(raw_data, quiet = TRUE)
      loaded <- TRUE
      break
    }
  }

  if (!loaded) {
    skip("No data files found for testing")
  }

  # Extract loci
  loci <- extract_loci(std_data, quiet = TRUE)

  # Check we found some loci
  expect_gt(length(loci), 0)

  # Check loci are valid HLA genes
  common_loci <- c("A", "B", "C", "DRB1", "DQB1", "DPB1")
  expect_true(any(loci %in% common_loci))
})

test_that("trim_hla_results works with real data", {
  # Load and standardize a data file
  formats <- c("csv", "tsv", "txt", "xlsx")
  loaded <- FALSE

  for (fmt in formats) {
    file_path <- file.path(data_dir, paste0("family_typing_data.", fmt))
    if (file.exists(file_path)) {
      raw_data <- load_typing_data(file_path, quiet = TRUE)
      std_data <- standardize_colnames(raw_data, quiet = TRUE)
      loaded <- TRUE
      break
    }
  }

  if (!loaded) {
    skip("No data files found for testing")
  }

  # Get allele columns
  allele_cols <- grep("_[12]$", names(std_data), value = TRUE)
  if (length(allele_cols) == 0) {
    skip("No valid allele columns found in data")
  }

  # Select one column with 4-field alleles for testing
  has_4fields <- FALSE
  for (col in allele_cols) {
    if (any(grepl("\\*[0-9]+:[0-9]+:[0-9]+:[0-9]+", std_data[[col]]))) {
      has_4fields <- TRUE
      break
    }
  }

  if (!has_4fields) {
    # If no 4-field alleles found, just check the function runs
    result_2field <- trim_hla_results(std_data, resolution = 2, quiet = TRUE)
    expect_is(result_2field, "data.frame")
  } else {
    # Test trimming to different resolutions
    result_2field <- trim_hla_results(std_data, resolution = 2, quiet = TRUE)

    # Check some columns were trimmed
    trimmed_found <- FALSE
    for (col in allele_cols) {
      orig_alleles <- std_data[[col]]
      trimmed_alleles <- result_2field[[col]]

      # Skip empty columns
      if (all(is.na(orig_alleles)) || all(orig_alleles == "")) next

      # Check if any allele was shortened
      if (any(nchar(trimmed_alleles, keepNA = TRUE) < nchar(orig_alleles, keepNA = TRUE))) {
        trimmed_found <- TRUE
        break
      }
    }

    expect_true(trimmed_found)
  }
})
