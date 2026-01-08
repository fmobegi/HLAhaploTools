# ~~~~~~~~~~~~~~~~~~~~~~~~
# calculate_hla_frequency
# ~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate HLA Allele Frequencies
#'
#' Computes allele frequencies for all HLA genes in a tidy dataset with `_1` and `_2` allele columns.
#' If no population column exists, a default "POP" value is used.
#' Handles multi-allele strings (e.g., "A*01:01:01:01/A*01:01:01:02") by splitting before counting.
#'
#' @param hped A data frame or tibble with HLA allele columns.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A tibble with columns: gene, allele, count, freq; sorted by gene and descending freq.
#' @examples
#' \dontrun{
#' # Calculate frequencies from formatted data
#' data <- load_typing_data("hla_data.csv") %>%
#'   reformat_typing_data()
#'
#' # Get allele frequencies
#' freqs <- calculate_hla_frequency(data)
#'
#' # Plot the frequencies
#' plot_hla_allele_frequency(data)
#' }
#'
#' @importFrom dplyr arrange desc bind_rows
#' @importFrom tibble tibble as_tibble
#'
#' @export
calculate_hla_frequency <- function(hped, quiet = FALSE) {
  if (!quiet) message("\t📊 Started Calculating HLA allele frequencies...")

  if (!is.data.frame(hped)) {
    stop("\t❌ Input must be a data frame or tibble.")
  }

  if (nrow(hped) == 0) {
    warning("\t⚠️ Input has zero rows. Returning empty frequency table.")
    tibble::tibble(
      gene = character(),
      allele = character(),
      count = integer(),
      freq = numeric(),
      pop = character()
    )
  }

  if (!"Population" %in% colnames(hped)) {
    if (!quiet) {
      message("\tℹ️ No population column found — assigning default 'POP'.")
    }
    hped$Population <- "POP"
  }

  allele_cols <- grep("_[12]$", names(hped), value = TRUE)
  if (length(allele_cols) == 0) {
    warning("\t❌ No HLA allele columns found (_1, _2 suffix).")
    tibble::tibble(
      gene = character(),
      allele = character(),
      count = integer(),
      freq = numeric(),
      pop = character()
    )
  }

  gene_names <- unique(sub("_[12]$", "", allele_cols))
  if (!quiet) {
    message(
      sprintf(
        "\t🧬 Found %d genes to analyze: %s",
        length(gene_names),
        paste(
          gene_names,
          collapse = ", "
        )
      )
    )
  }

  process_one_gene <- function(gene) {
    tryCatch(
      {
        gene_cols <- paste0(gene, c("_1", "_2"))
        gene_cols <- gene_cols[gene_cols %in% names(hped)]

        alleles <- unlist(hped[gene_cols],
          use.names = FALSE
        )
        alleles <- alleles[!is.na(alleles) & alleles != ""]
        alleles <- trimws(alleles)

        # Exclude ambiguous multi-allele calls (those with "/")
        alleles_clean <- alleles[!grepl("/", alleles)]

        if (length(alleles_clean) == 0) {
          if (!quiet) message(sprintf("\t⚠️ Gene '%s' has no valid allele data — skipped.", gene))
          NULL
        }

        counts <- table(alleles_clean)
        result <- tibble::tibble(
          gene = gene,
          allele = names(counts),
          count = as.integer(counts),
          freq = as.integer(counts) / length(alleles_clean),
          pop = unique(hped$Population)[1]
        )

        if (!quiet) {
          message(sprintf(
            "\t✅  %s has %d unique alleles from %d total.",
            gene, length(counts), length(alleles_clean)
          ))
        }

        result
      },
      error = function(e) {
        warning(sprintf("\t❌ Skipping '%s' due to error: %s", gene, e$message))
        NULL
      }
    )
  }

  results_list <- lapply(gene_names, process_one_gene)
  df_freq <- dplyr::bind_rows(results_list)

  if (nrow(df_freq) == 0) {
    warning("\t⚠️ No valid allele data found.")
    tibble::tibble(
      gene = character(),
      allele = character(),
      count = integer(),
      freq = numeric(),
      pop = character()
    )
  }

  df_freq <- dplyr::arrange(df_freq, gene, dplyr::desc(freq))

  if (!quiet) {
    message(sprintf(
      "\t🏁 Completed frequency analysis: %d alleles across %d loci.",
      nrow(df_freq), length(unique(df_freq$gene))
    ))
  }

  df_freq
}
