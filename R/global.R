#' Calculate HLA Allele Frequencies
#'
#' Computes allele frequencies for all HLA genes in a tidy dataset with `_1` and `_2` allele columns.
#' If no population column exists, a default "FSH_Immuno" value is used.
#' Handles multi-allele strings (e.g., "A*01:01:01:01/A*01:01:01:02") by splitting before counting.
#'
#' @param hped A data frame or tibble with HLA allele columns.
#'
#' @return A tibble with columns: gene, allele, count, freq; sorted by gene and descending freq.
#' @importFrom dplyr arrange desc
#' @importFrom tibble tibble as_tibble
#' @export
calculate_HLA_frequency <- function(hped) {
  message("\t📊 Calculating HLA allele frequencies...")

  if (!"Population" %in% colnames(hped)) {
    message("\tℹ️ No population column found — assigning default 'FSH_Immuno'.")
    hped$Population <- "FSH"
  }

  allele_cols <- grep("_[12]$", names(hped), value = TRUE)
  gene_names <- unique(sub("_[12]$", "", allele_cols))

  df_list <- lapply(gene_names, function(gene) {
    tryCatch(
      {
        gene_cols <- paste0(gene, c("_1", "_2"))
        gene_cols <- gene_cols[gene_cols %in% names(hped)]

        alleles <- unlist(hped[gene_cols], use.names = FALSE)
        alleles <- alleles[!is.na(alleles)]

        # Split multi-allele strings (e.g. A*01:01/A*01:02) into separate observations
        alleles_split <- unlist(strsplit(alleles, "/"))
        alleles_split <- trimws(alleles_split)

        if (length(alleles_split) == 0) {
          message(sprintf("\t⚠️ Gene '%s' has no allele data — skipped.", gene))
          return(NULL)
        }

        counts <- table(alleles_split)
        tibble::tibble(
          gene = gene,
          allele = names(counts),
          count = as.integer(counts),
          freq = as.integer(counts) / length(alleles_split)
        )
      },
      error = function(e) {
        warning(sprintf("\t❌ Skipping gene '%s' due to error: %s", gene, e$message))
        NULL
      }
    )
  })

  freq_result <- do.call(rbind, Filter(Negate(is.null), df_list))
  message("\t✅ Allele frequency calculation complete.\n")

  return(freq_result %>% dplyr::arrange(gene, dplyr::desc(freq)))
}


# ~~~~~~~~~~~~~~~~~~~ PLOT FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Plot HLA Allele Frequencies by Gene
#'
#' Visualizes allele frequencies across genes, grouping rare alleles into "others".
#' Frequencies are computed using `calculate_HLA_frequency()`.
#'
#' @param hped A data frame or tibble with HLA allele columns.
#' @param minFreq Minimum frequency threshold (default = 0.05) above which allele labels are shown.
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_fill_manual xlab ylab theme_classic theme
#' @importFrom dplyr mutate if_else group_by summarise desc
#' @importFrom tibble as_tibble
#' @importFrom RColorBrewer brewer.pal
#' @export
Plot_HLA_allele_frequency <- function(hped, minFreq = 0.05) {
  freq_file <- tryCatch(calculate_HLA_frequency(hped), error = function(e) {
    stop("Failed to calculate frequencies: ", e$message)
  })

  freq_plot_df <- as_tibble(freq_file) %>%
    dplyr::mutate(Label = dplyr::if_else(freq > minFreq, allele, "others"))

  custom_colors <- c(
    colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(freq_plot_df$Label))),
    "gray90"
  )

  freq_plot_df %>%
    dplyr::group_by(gene, Label) %>%
    dplyr::summarise(freq_sum = sum(freq), .groups = "drop") %>%
    ggplot2::ggplot(aes(x = gene, y = freq_sum, fill = Label)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(aes(label = Label), position = position_stack(0.5), size = 2.5) +
    ggplot2::theme_classic() +
    ggplot2::xlab("HLA Genes") +
    ggplot2::ylab("Allele Frequency") +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::theme(legend.position = "none")
} #' Plot HLA Allele Frequencies by Gene
#'
#' Visualizes allele frequencies across genes, grouping rare alleles into "others".
#' Frequencies are computed using `calculate_HLA_frequency()`.
#'
#' @param hped A data frame or tibble with HLA allele columns.
#' @param minFreq Minimum frequency threshold (default = 0.05) above which allele labels are shown.
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_fill_manual xlab ylab theme_classic theme
#' @importFrom dplyr mutate if_else group_by summarise desc
#' @importFrom tibble as_tibble
#' @importFrom RColorBrewer brewer.pal
#' @export
Plot_HLA_allele_frequency <- function(hped, minFreq = 0.05) {
  freq_file <- tryCatch(calculate_HLA_frequency(hped), error = function(e) {
    stop("Failed to calculate frequencies: ", e$message)
  })

  freq_plot_df <- as_tibble(freq_file) %>%
    dplyr::mutate(Label = dplyr::if_else(freq > minFreq, allele, "others"))

  custom_colors <- c(
    colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(freq_plot_df$Label))),
    "gray90"
  )

  freq_plot_df %>%
    dplyr::group_by(gene, Label) %>%
    dplyr::summarise(freq_sum = sum(freq), .groups = "drop") %>%
    ggplot2::ggplot(aes(x = gene, y = freq_sum, fill = Label)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(aes(label = Label), position = position_stack(0.5), size = 2.5) +
    ggplot2::theme_classic() +
    ggplot2::xlab("HLA Genes") +
    ggplot2::ylab("Allele Frequency") +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::theme(legend.position = "none")
}


# ----------------------------------------------------------------------------------

#' Plot Count of Unique Alleles per HLA Gene
#'
#' Visualizes the number of unique observed alleles for each gene.
#'
#' @param hped A data frame or tibble with HLA allele columns.
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_classic scale_fill_manual xlab ylab theme element_text
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr group_by summarise
#' @export
Plot_HLA_allele_count <- function(hped) {
  freq_file <- tryCatch(calculate_HLA_frequency(hped), error = function(e) {
    stop("Failed to calculate frequencies: ", e$message)
  })

  count_df <- freq_file %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(Allele_Count = dplyr::n(), .groups = "drop")

  ggplot2::ggplot(count_df, aes(x = gene, y = Allele_Count, fill = gene)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::xlab("HLA Gene") +
    ggplot2::ylab("Unique Alleles Observed") +
    ggplot2::scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(nrow(count_df))) +
    ggplot2::theme(
      legend.position = "none",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16)
    )
}

# ---------------------------------------------------------------------------


#' Plot HLA Allele Diversity by Population
#'
#' Shows the relative frequency of top alleles for a specific gene within each population.
#'
#' @param hped A data frame or tibble with paired `_1` / `_2` allele columns and a population identifier.
#' @param gene The HLA gene to analyze (e.g., "A", "DRB1").
#' @param ntop Number of top alleles to label per population.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate summarise group_by arrange filter if_else
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_fill_manual theme_classic theme labs
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_HLA_Diversity <- function(hped, gene = "A", ntop = 10) {
  message("\t📊 Preparing HLA diversity plot for gene: ", gene)

  # Ensure Population column is present
  if (!"Population" %in% names(hped)) {
    message("\tℹ️ No 'Population' column found — assigning all records to 'FSH_Immuno'")
    hped$Population <- "FSH"
  }

  # Reshape to long format
  allele_cols <- grep("_[12]$", names(hped), value = TRUE)
  long_df <- hped %>%
    tidyr::pivot_longer(cols = all_of(allele_cols), names_to = "Column", values_to = "allele") %>%
    dplyr::mutate(HLA_gene = sub("_.*", "", Column)) %>%
    dplyr::filter(!is.na(allele)) %>%
    # Split MAC multi-allele strings
    tidyr::separate_rows(allele, sep = "/")

  # Group and summarize frequencies
  allele_summary <- long_df %>%
    dplyr::group_by(Population, HLA_gene, allele) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(Population, HLA_gene) %>%
    dplyr::mutate(
      freq = count / sum(count),
      Label = dplyr::if_else(freq > 0.05, allele, "others")
    )

  # Identify top alleles
  top_alleles <- allele_summary %>%
    dplyr::filter(HLA_gene == gene) %>%
    dplyr::arrange(Population, dplyr::desc(freq)) %>%
    dplyr::group_by(Population) %>%
    dplyr::slice_head(n = ntop) %>%
    dplyr::pull(allele) %>%
    unique()

  # Prep for plotting
  plot_df <- allele_summary %>%
    dplyr::filter(HLA_gene == gene) %>%
    dplyr::mutate(Label = dplyr::if_else(allele %in% top_alleles, allele, "others")) %>%
    dplyr::group_by(Population, Label) %>%
    dplyr::summarise(freq2 = sum(freq), .groups = "drop")

  ggplot2::ggplot(plot_df, aes(x = Population, y = freq2, fill = Label)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4.5) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Population", y = "Allele Frequency") +
    ggplot2::scale_fill_manual(values = c(
      colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(plot_df$Label))),
      "gray90"
    )) +
    ggplot2::theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.position = "none"
    )
}

## PAckages ----
# usethis::use_package("janitor")
# usethis::use_package("dplyr")
# usethis::use_package("tidyr")
# usethis::use_package("purrr")
# usethis::use_package("tibble")
# usethis::use_package("forcats")
# usethis::use_package("readr")
# usethis::use_package("readxl")
# usethis::use_package("future")
# usethis::use_package("furrr")
# usethis::use_package("future.apply")
# usethis::use_package("immunotation")
# usethis::use_package("HLAtools")
