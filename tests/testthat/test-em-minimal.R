# Test EM algorithm with minimal synthetic data

test_that("EM algorithm works with minimal synthetic data", {
  # Load the minimal test data
  test_file <- "testdata/mini_test_data.csv"
  if (!file.exists(test_file)) {
    skip("Minimal test data file not found")
  }

  # Load and prepare data
  test_data <- load_typing_data(test_file, quiet = TRUE)
  expect_equal(nrow(test_data), 5) # Verify we have 5 samples

  # Standardize column names
  std_data <- standardize_colnames(test_data, quiet = TRUE)

  # Extract loci
  loci <- extract_loci(std_data, quiet = TRUE)
  expect_equal(sort(loci), c("A", "B"))

  # Collapse genotypes
  genotypes <- collapse_genotypes(std_data, loci, quiet = TRUE)
  expect_equal(nrow(genotypes), 5)

  # Enumerate all possible diplotypes
  all_diplotypes <- list()
  for (i in seq_len(nrow(genotypes))) {
    diplotypes <- enumerate_diplotypes(genotypes[i, , drop = FALSE], quiet = TRUE)
    all_diplotypes[[i]] <- diplotypes
  }

  # Run EM algorithm
  hap_freqs <- run_em_algorithm(all_diplotypes, max_iter = 10, quiet = TRUE)

  # Check results
  expect_true(is.numeric(hap_freqs))
  expect_true(all(hap_freqs >= 0))
  expect_true(all(hap_freqs <= 1))
  expect_equal(sum(hap_freqs), 1, tolerance = 1e-6)

  # Expected haplotypes should include these combinations
  expected_haplotypes <- c(
    "A*01:01~B*07:02", "A*01:01~B*08:01",
    "A*02:01~B*07:02", "A*02:01~B*08:01",
    "A*02:01~B*44:02", "A*03:01~B*44:02"
  )

  # Check that at least some expected haplotypes are found
  found_expected <- sum(expected_haplotypes %in% names(hap_freqs))
  expect_gt(found_expected, 0)

  # Run posterior computation
  posteriors <- compute_posteriors(all_diplotypes, hap_freqs, quiet = TRUE)
  expect_equal(length(posteriors), length(all_diplotypes))

  # Check that posteriors sum to 1 for each individual
  for (i in seq_along(posteriors)) {
    if (length(posteriors[[i]]) > 0) {
      expect_equal(sum(posteriors[[i]]), 1, tolerance = 1e-6)
    }
  }
})

test_that("infer_haplotypes works with minimal synthetic data", {
  # Load the minimal test data
  test_file <- "testdata/mini_test_data.csv"
  if (!file.exists(test_file)) {
    skip("Minimal test data file not found")
  }

  # Run the full inference pipeline
  results <- infer_haplotypes(
    df = read.csv(test_file),
    parallel = FALSE,
    quiet = TRUE,
    isfamily = FALSE
  )

  # Check results structure
  expect_true(is.list(results))
  expect_true("haplotype_frequencies" %in% names(results))
  expect_true("top_diplotypes" %in% names(results))
  expect_true("posteriors" %in% names(results))
  expect_true("loci_used" %in% names(results))

  # Check that frequencies sum to 1
  hap_freqs <- results$haplotype_frequencies$frequency
  expect_equal(sum(hap_freqs), 1, tolerance = 1e-6)

  # Check that all samples have a top diplotype assigned
  expect_equal(nrow(results$top_diplotypes), 5)

  # Check that loci were correctly identified
  expect_equal(sort(results$loci_used), c("A", "B"))
})
