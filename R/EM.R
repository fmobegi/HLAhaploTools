#' Standardize HLA Genotype Column Names
#'
#' Converts genotype column names (e.g. A.1) to A_1 format and applies consistent styling.
#'
#' @param df A raw HLA genotype data frame.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A cleaned data frame with standardized column names.
#' @examples
#' \dontrun{
#' raw_data <- read.csv("hla_data.csv")
#' clean_data <- standardize_colnames(raw_data)
#' }
#' @export
#' @importFrom janitor clean_names
standardize_colnames <- function(df, quiet = FALSE) {
  if (!quiet) message("\t🔧 Standardizing genotype column names...")
  if (!is.data.frame(df)) stop("\t❌ Input must be a data frame or tibble.")

  names(df) <- gsub("[\\.\\- ]", "_", names(df))
  df <- janitor::clean_names(df, case = "none")

  if (!quiet) {
    hla_cols <- sum(grepl("_[12]$", names(df)))
    message(sprintf(
      "\t✅ Standardized %d column names (%d HLA allele columns).",
      length(names(df)), hla_cols
    ))
  }
  df
}

#' Extract HLA Locus Names from Genotype Table
#'
#' Identifies unique locus names based on "_1" and "_2" suffixes.
#'
#' @param df A data frame with standardized column names.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A character vector of distinct locus names (e.g. A, B, DRB1).
#' @examples
#' \dontrun{
#' loci <- extract_loci(hla_data)
#' print(loci)
#' }
#' @export
extract_loci <- function(df, quiet = FALSE) {
  if (!quiet) message("\t🔬 Extracting HLA loci from column names...")
  if (!is.data.frame(df)) stop("\t❌ Input must be a data frame or tibble.")

  allele_cols <- grep("_[12]$", names(df), value = TRUE)
  if (length(allele_cols) == 0) {
    warning("\t⚠️ No HLA allele columns found with _1 or _2 suffix.")
    return(character(0))
  }

  loci <- unique(gsub("_[12]$", "", allele_cols))
  if (!quiet) {
    message(sprintf(
      "\t✅ Found %d HLA loci: %s",
      length(loci), paste(loci, collapse = ", ")
    ))
  }
  loci
}

#' Run haplotype EM estimation and summarise results
#'
#' Reads a family typing data file, runs the EM algorithm using
#' \code{haplo.stats::haplo.em}, and returns a tidy tibble of haplotypes.
#'
#' @param file Path to the family typing data file (TSV format).
#' @param loci Character vector of locus names to include.
#' @param sep Separator string used when concatenating alleles into haplotype strings.
#'   Defaults to \code{"~"}.
#'
#' @return A tibble with three columns: Haplotype, EM_Probability, Frequency.
#' @examples
#' \dontrun{
#' loci <- c("A", "B", "C", "DRB1", "DQB1", "DPB1", "DQA1", "DPA1", "DRB4", "G", "H", "MICA", "MICB")
#' hap_summary <- run_haplo_em("./inst/extdata/family_typing_data.tsv", loci)
#' print(hap_summary)
#' }
#' @import dplyr
#' @import tibble
#' @import haplo.stats
#' @importFrom stats na.omit
#' @export
run_haplo_em <- function(file, loci, sep = "~") {
  cols <- c(paste0(loci, "_1"), paste0(loci, "_2"))

  df <- utils::read.delim(file) %>%
    dplyr::select(FAMILY_ID, dplyr::all_of(cols)) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(-FAMILY_ID, ~ dplyr::na_if(.x, ""))) %>%
    dplyr::select(-c(FAMILY_ID)) %>%
    as.matrix()

  em <- haplo.stats::haplo.em(df, locus.label = loci, miss.val = NA)

  hap_df <- data.frame(em$haplotype, em$hap.prob) %>%
    tibble::as_tibble() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Haplotype = paste(stats::na.omit(dplyr::c_across(1:(ncol(.) - 1))), collapse = sep)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Haplotype, EM_Probability = em.hap.prob) %>%
    dplyr::group_by(Haplotype) %>%
    dplyr::summarise(
      EM_Probability = sum(EM_Probability),
      Frequency = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Haplotype, dplyr::desc(EM_Probability))

  hap_df
}

#' Plot Top HLA Haplotypes by Frequency
#'
#' Visualizes the most common haplotypes and their frequencies.
#'
#' @param hap_freq A tibble with haplotype frequencies.
#' @param n_top Number of top haplotypes to display (default = 20).
#' @param min_freq Minimum frequency to include (default = 0.01).
#' @param color Fill color for bars (default = "steelblue").
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' result <- run_haplo_em("family_typing_data.tsv", loci)
#' p <- plot_top_haplotypes(result, n_top = 10, min_freq = 0.02)
#' }
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_y_continuous theme_minimal labs element_text coord_flip theme ggtitle theme_void
#' @export
plot_top_haplotypes <- function(hap_freq,
                                n_top = 20,
                                min_freq = 0.01,
                                color = "steelblue",
                                quiet = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("\t❌ Package ggplot2 is required. Install with install.packages('ggplot2')")
  }

  if (!quiet) message("\t📊 Creating haplotype frequency plot...")

  if (!all(c("Haplotype", "Frequency") %in% names(hap_freq))) {
    stop("\t❌ Input must have 'Haplotype' and 'Frequency' columns.")
  }

  plot_data <- hap_freq[hap_freq$Frequency >= min_freq, ]
  plot_data <- plot_data[order(plot_data$Frequency, decreasing = TRUE), ]

  if (nrow(plot_data) == 0) {
    warning(sprintf("\t⚠️ No haplotypes meet minimum frequency threshold of %.4f", min_freq))
    return(ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::ggtitle("No haplotypes meet frequency threshold"))
  }

  if (nrow(plot_data) > n_top) {
    plot_data <- plot_data[seq_len(n_top), ]
    if (!quiet) {
      message(sprintf("\t✓ Showing top %d haplotypes.", n_top))
    }
  }

  plot_data$Haplotype <- factor(plot_data$Haplotype, levels = rev(plot_data$Haplotype))

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Haplotype, y = Frequency)) +
    ggplot2::geom_bar(stat = "identity", fill = color, width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", Frequency * 100)), hjust = -0.2) +
    ggplot2::scale_y_continuous(
      labels = scales::percent,
      limits = c(0, max(plot_data$Frequency) * 1.2)
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Top HLA Haplotype Frequencies", x = "Haplotype", y = "Frequency") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}
