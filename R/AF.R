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
    return(tibble::tibble(
      gene = character(),
      allele = character(),
      count = integer(),
      freq = numeric(),
      pop = character()
    ))
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
    return(tibble::tibble(
      gene = character(),
      allele = character(),
      count = integer(),
      freq = numeric(),
      pop = character()
    ))
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
          return(NULL)
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

        return(result)
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
    return(tibble::tibble(
      gene = character(),
      allele = character(),
      count = integer(),
      freq = numeric(),
      pop = character()
    ))
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~
# plot_hla_allele_frequency
# ~~~~~~~~~~~~~~~~~~~~~~~~~
#' Plot HLA Allele Frequencies by Gene
#'
#' Visualizes allele frequencies across genes, grouping rare alleles into "others".
#' Automatically splits classical and non-classical HLA genes into separate panels with nested facets.
#'
#' @param freq_data A data frame with gene, allele, and freq columns.
#' @param min_freq Minimum frequency threshold (default = 0.05) above which allele labels are shown.
#' @param split_threshold Threshold for number of unique genes that triggers multi-panel layout (default = 11).
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A ggplot object that can be printed or composed with patchwork.
#'
#' @import ggplot2
#' @importFrom dplyr group_by mutate summarise arrange desc filter case_when
#' @importFrom tidyr replace_na
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom ggh4x facet_nested
#'
#' @export
plot_hla_allele_frequency <- function(freq_data,
                                      min_freq = 0.05,
                                      split_threshold = 11,
                                      quiet = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("❌ Please install ggplot2")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("❌ Please install patchwork")
  }
  if (!requireNamespace("ggh4x", quietly = TRUE)) {
    stop("❌ Please install ggh4x")
  }

  if (!quiet) {
    message("\t📊 Starting allele frequency plot generation...")
  }

  req_cols <- c("gene", "allele", "freq")
  if (!all(req_cols %in% colnames(freq_data))) {
    stop(
      "\t❌ Missing required columns: ",
      paste(setdiff(
        req_cols,
        colnames(freq_data)
      ), collapse = ", ")
    )
  }

  # Gene classification
  class_i <- c("A", "B", "C")
  class_ii <- c("DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5")
  nonclass_i <- c("E", "F", "G", "H", "J", "MICA", "MICB", "HFE")
  nonclass_ii <- c("DMA", "DMB", "DOA", "DOB")

  freq_data <- freq_data %>%
    dplyr::mutate(
      super_class = dplyr::case_when(
        gene %in% c(class_i, class_ii) ~ "Classical",
        gene %in% c(nonclass_i, nonclass_ii) ~ "Nonclassical",
        TRUE ~ "Other"
      ),
      labels = dplyr::case_when(
        gene %in% class_i ~ "Class I",
        gene %in% class_ii ~ "Class II",
        gene %in% nonclass_i ~ "Class I",
        gene %in% nonclass_ii ~ "Class II",
        TRUE ~ "Other"
      )
    )

  unique_genes <- unique(freq_data$gene)
  need_split <- length(unique_genes) > split_threshold

  # Show genes present by class
  classical_genes_present <- sort(
    unique(freq_data$gene[freq_data$super_class == "Classical"])
  )
  nonclassical_genes_present <- sort(
    unique(freq_data$gene[freq_data$super_class == "Nonclassical"])
  )

  if (!quiet) {
    message(
      "\t🔍 Classical genes present: ",
      paste(classical_genes_present, collapse = ", ")
    )
    message(
      "\t🔍 Nonclassical genes present: ",
      paste(nonclassical_genes_present, collapse = ", ")
    )
    message(
      sprintf(
        "\t🔢 Total unique genes found: %d (split threshold: %d)",
        length(unique_genes),
        split_threshold
      )
    )
    if (need_split) {
      message("\t⚙️ Plot will be split into Classical and Nonclassical panels.")
    } else {
      message("\t⚙️ Plotting all genes in a single panel.")
    }
  }

  # Faceted panel builder
  create_nested_panel <- function(data, title = NULL) {
    if (nrow(data) == 0) {
      ggplot2::ggplot() +
        ggplot2::theme_void()
    }

    if (!quiet) {
      message(sprintf(
        "\t\t🎨 Creating panel for %d rows of data...",
        nrow(data)
      ))
    }

    plot_data <- data %>%
      dplyr::group_by(gene) %>%
      dplyr::mutate(
        allele_label = ifelse(freq >= min_freq, allele, "Other"),
        allele_label = tidyr::replace_na(allele_label, "Other")
      ) %>%
      dplyr::group_by(super_class, labels, gene, allele_label) %>%
      dplyr::summarise(freq = sum(freq), .groups = "drop") %>%
      dplyr::arrange(gene, dplyr::desc(freq))

    gene_colors <- create_gene_palette(unique(plot_data$gene))
    color_map <- create_freq_gradient_palette(plot_data, gene_colors)

    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = gene,
        y = freq,
        fill = interaction(gene, allele_label)
      )
    ) +
      ggplot2::geom_bar(
        stat = "identity",
        position = "stack",
        color = "white",
        linewidth = 0.2
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = ifelse(
          freq >= min_freq * 2,
          sprintf("%s\n(%.1f%%)", allele_label, freq * 100),
          ""
        )),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 2.5
      ) +
      ggplot2::scale_fill_manual(values = color_map) +
      ggh4x::facet_nested(. ~ super_class + labels,
        scales = "free",
        space = "free"
      ) +
      ggplot2::labs(title = title, y = "Frequency", x = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = 12
        ),
        strip.text = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
  }

  if (need_split) {
    if (!quiet) {
      message("\t📑 Creating nested panels for Classical and Nonclassical genes...")
    }

    classical_data <- dplyr::filter(
      freq_data,
      super_class == "Classical"
    )
    nonclassical_data <- dplyr::filter(
      freq_data,
      super_class == "Nonclassical"
    )

    panels <- list()
    if (nrow(classical_data) > 0) {
      panels <- c(panels, list(create_nested_panel(classical_data)))
    }
    if (nrow(nonclassical_data) > 0) {
      panels <- c(panels, list(create_nested_panel(nonclassical_data)))
    }

    if (length(panels) == 2) {
      p <- panels[[1]] / panels[[2]] +
        patchwork::plot_layout(ncol = 1, heights = c(3, 2))
    } else {
      p <- panels[[1]]
    }

    p <- p + patchwork::plot_annotation(
      title = "HLA Allele Frequencies by Gene Class",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          size = 14,
          face = "bold"
        )
      )
    )
  } else {
    if (!quiet) {
      message("\t🖼️ Plotting all genes in one panel")
    }
    p <- create_nested_panel(freq_data, "HLA Allele Frequencies by Gene")
  }

  if (!quiet) {
    message("\t✅ Allele Frequency Plot created successfully.")
  }
  return(p)
}

# ~~~~~~~~~~~~~~~~~~~~~
# plot_hla_allele_count
# ~~~~~~~~~~~~~~~~~~~~~
#' Plot Count of Unique Alleles per HLA Gene
#'
#' Visualizes the number of unique observed alleles for each gene.
#'
#' @param freq_data A data frame with gene and allele columns.
#' @param min_count Minimum unique allele count to include gene (default = 0).
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A ggplot object that can be further customized or printed.
#'
#' @import ggplot2
#' @import dplyr

#' @export
plot_hla_allele_count <- function(freq_data,
                                  min_count = 0,
                                  quiet = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "\t❌ Package ggplot2 is required for plotting. Please install it with install.packages('ggplot2')"
    )
  }

  if (!quiet) {
    message("\t📊 Starting allele count plot generation...")
  }

  # Validate input data
  req_cols <- c("gene", "allele")
  if (!all(req_cols %in% colnames(freq_data))) {
    stop("\t❌ Input data must contain columns: gene, allele")
  }

  # Gene classification to report classes found
  class_i <- c("A", "B", "C")
  class_ii <- c("DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5")
  nonclass_i <- c("E", "F", "G", "H", "J", "MICA", "MICB", "HFE")
  nonclass_ii <- c("DMA", "DMB", "DOA", "DOB")

  # Count unique alleles per gene
  if (!quiet) {
    message("\t🔍 Counting unique alleles per gene...")
  }
  allele_counts <- freq_data %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      unique_alleles = dplyr::n_distinct(allele),
      .groups = "drop"
    ) %>%
    dplyr::filter(unique_alleles >= min_count) %>%
    dplyr::arrange(dplyr::desc(unique_alleles))

  if (nrow(allele_counts) == 0) {
    warning(sprintf("\t⚠️ No genes with at least %d unique alleles.", min_count))
    return(ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::ggtitle(sprintf(
        "No genes with ≥%d unique alleles", min_count
      )))
  }

  # Determine classes found in data
  found_genes <- allele_counts$gene
  classical_genes <- sort(unique(
    found_genes[found_genes %in%
      c(class_i, class_ii)]
  ))
  nonclassical_genes <- sort(unique(
    found_genes[found_genes %in%
      c(nonclass_i, nonclass_ii)]
  ))
  other_genes <- sort(unique(
    found_genes[!found_genes %in%
      c(class_i, class_ii, nonclass_i, nonclass_ii)]
  ))

  if (!quiet) {
    message("\t🧬 Genes found in each class:")
    message("\t\tClassical genes: ", ifelse(
      length(classical_genes) > 0,
      paste(classical_genes, collapse = ", "),
      "None"
    ))
    message(
      "\t\tNonclassical genes: ",
      ifelse(
        length(nonclassical_genes) > 0,
        paste(nonclassical_genes, collapse = ", "),
        "None"
      )
    )
    message("\t\tOther genes: ", ifelse(
      length(other_genes) > 0,
      paste(other_genes, collapse = ", "),
      "None"
    ))
  }

  if (!quiet) {
    message(sprintf(
      "\t🎨 Preparing plot for %d genes...",
      nrow(allele_counts)
    ))
  }
  gene_colors <- create_gene_palette(allele_counts$gene)

  p <- ggplot2::ggplot(
    allele_counts,
    ggplot2::aes(
      x = reorder(gene, unique_alleles),
      y = unique_alleles,
      fill = gene
    )
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = unique_alleles),
      position = ggplot2::position_stack(vjust = 0.5),
      color = "white"
    ) +
    ggplot2::scale_fill_manual(values = gene_colors) +
    ggplot2::ylab("Number of Unique Alleles") +
    ggplot2::xlab("HLA Gene") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(
        face = "bold",
        angle = 45,
        hjust = 1
      ),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::ggtitle("HLA Allelic Diversity")

  if (!quiet) {
    message("\t✅ Allele count plot created successfully.")
  }
  return(p)
}

# ~~~~~~~~~~~~~~~~~~
# plot_hla_diversity
# ~~~~~~~~~~~~~~~~~~
#' Plot HLA Allele Diversity by Population
#'
#' Shows the relative frequency of top alleles for a specific gene within each population.
#'
#' @param freq_data A data frame with gene, allele, freq and pop columns.
#' @param gene The HLA gene to analyze (e.g., "A", "DRB1").
#' @param ntop Number of top alleles to label per population (default = 5).
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A ggplot object that can be further customized or printed.
#'
#' @import ggplot2
#' @importFrom forcats fct_lump_n
#' @importFrom dplyr group_by mutate summarise arrange desc
#'
#' @export
plot_hla_diversity <- function(freq_data,
                               gene = "A",
                               ntop = 5,
                               quiet = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "\t❌ Package 'ggplot2' is required. Please install it with install.packages('ggplot2')"
    )
  }
  if (!requireNamespace("forcats", quietly = TRUE)) {
    stop(
      "\t❌ Package 'forcats' is required. Please install it with install.packages('forcats')"
    )
  }

  if (!quiet) {
    message(sprintf("\t📊 Generating allele diversity plot for gene %s...", gene))
  }

  # Validate input data
  req_cols <- c("gene", "allele", "freq", "pop")
  if (!all(req_cols %in% colnames(freq_data))) {
    stop("\t❌ Input data must contain columns: gene, allele, freq, pop")
  }

  # Filter data for the selected gene
  gene_data <- freq_data[freq_data$gene == gene, ]

  if (nrow(gene_data) == 0) {
    warning(sprintf("\t⚠️ No allele data found for gene '%s'.", gene))
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::ggtitle(sprintf("No allele data for %s", gene))
    )
  }

  if (!quiet) {
    message("\t🔍 Filtering top ", ntop, " alleles per population...")
  }

  # For each population, keep only top N alleles and group the rest as "Other"
  plot_data <- gene_data %>%
    dplyr::group_by(pop) %>%
    dplyr::mutate(allele_group = forcats::fct_lump_n(
      allele,
      n = ntop,
      w = freq,
      other_level = "Other"
    )) %>%
    dplyr::group_by(pop, allele_group) %>%
    dplyr::summarise(freq = sum(freq), .groups = "drop") %>%
    dplyr::arrange(pop, dplyr::desc(freq))

  # Generate color palette specific to population+allele combinations
  color_map <- create_pop_allele_palette(plot_data, gene)

  # Compose the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = pop,
    y = freq,
    fill = interaction(pop, allele_group)
  )) +
    ggplot2::geom_bar(
      stat = "identity",
      position = "stack",
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(
        freq >= 0.05,
        sprintf("%s\n(%.1f%%)", allele_group, freq * 100),
        ""
      )),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 2.8
    ) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::labs(
      title = sprintf("%s Allele Distribution by Population", gene),
      x = "Population",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(
        face = "bold",
        angle = 45,
        hjust = 1
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 13,
        face = "bold"
      )
    )

  if (!quiet) {
    message("\t✅ Diversity plot created successfully.")
  }
  return(p)
}

## --------------- END OF EXPORTED FUNCTIONS ------------------------

# Helper function to create a palette of distinct colors for genes
create_gene_palette <- function(genes) {
  # Predefined colors for common HLA genes
  hla_colors <- c(
    "A" = "#E41A1C", # Red
    "B" = "#377EB8", # Blue
    "C" = "#4DAF4A", # Green
    "DRB1" = "#984EA3", # Purple
    "DRB3" = "#925E93", # Light Purple
    "DRB4" = "#7F4E83", # Medium Purple
    "DRB5" = "#6C3E73", # Dark Purple
    "DQB1" = "#FF7F00", # Orange
    "DQA1" = "#FFFF33", # Yellow
    "DPB1" = "#A65628", # Brown
    "DPA1" = "#F781BF", # Pink
    "E" = "#999999", # Gray
    "F" = "#66C2A5", # Teal
    "G" = "#FC8D62", # Salmon
    "H" = "#8DA0CB", # Light Blue
    "J" = "#E78AC3", # Orchid
    "MICA" = "#A6D854", # Light Green
    "MICB" = "#FFD92F", # Gold
    "HFE" = "#B3B3B3", # Silver
    "DMA" = "#B2DF8A", # Pale Green
    "DMB" = "#FB9A99", # Pale Red
    "DOA" = "#CAB2D6", # Lavender
    "DOB" = "#FDBF6F" # Light Orange
  )

  result <- character(length(genes))
  names(result) <- genes

  # Track how many dynamic colors have been generated
  dynamic_index <- 1

  for (i in seq_along(genes)) {
    gene <- genes[i]
    if (gene %in% names(hla_colors)) {
      result[i] <- hla_colors[gene]
    } else {
      # Use golden ratio-based HSV spacing for new colors
      h <- (dynamic_index * 0.618033988749895) %% 1
      result[i] <- hsv(h, 0.6, 0.9) # softer saturation and value
      dynamic_index <- dynamic_index + 1
    }
  }

  result
}

# Helper function to create gradient palettes based on frequency
create_freq_gradient_palette <- function(data, gene_colors) {
  result <- character()

  # For each gene-allele combination, create color
  for (gene in unique(data$gene)) {
    gene_data <- data[data$gene == gene, ]
    base_color <- gene_colors[gene]

    # Sort by frequency for consistent gradient
    gene_data <- gene_data[order(gene_data$freq, decreasing = TRUE), ]

    # Create gradient
    n_alleles <- nrow(gene_data)

    # Convert base color to RGB
    rgb_base <- col2rgb(base_color) / 255

    for (i in 1:n_alleles) {
      # Calculate fade factor (0 = original color, 0.8 = almost white)
      fade <- min(0.8, (i - 1) / max(1, n_alleles - 1))

      # Mix with white based on fade factor
      r <- rgb_base[1, 1] + (1 - rgb_base[1, 1]) * fade
      g <- rgb_base[2, 1] + (1 - rgb_base[2, 1]) * fade
      b <- rgb_base[3, 1] + (1 - rgb_base[3, 1]) * fade

      color <- rgb(r, g, b)

      # Create unique identifier for interaction term
      if ("allele_label" %in% colnames(gene_data)) {
        id <- paste(gene, gene_data$allele_label[i], sep = ".")
      } else if ("allele_group" %in% colnames(gene_data)) {
        id <- paste(gene_data$pop[i], gene_data$allele_group[i], sep = ".")
      } else {
        id <- gene
      }

      result[id] <- color
    }
  }

  result
}

# For plot_hla_diversity, modify the color creation
create_pop_allele_palette <- function(data, gene) {
  # Get a base color for this gene
  gene_colors <- create_gene_palette(gene)
  base_color <- gene_colors[[1]]
  result <- character()

  # For each population
  for (pop in unique(data$pop)) {
    pop_data <- data[data$pop == pop, ]

    # Sort by frequency
    pop_data <- pop_data[order(pop_data$freq, decreasing = TRUE), ]

    # Create gradient
    n_alleles <- nrow(pop_data)
    rgb_base <- col2rgb(base_color) / 255

    for (i in 1:n_alleles) {
      # More fading for lower frequencies
      fade <- min(0.8, (i - 1) / max(1, n_alleles - 1))

      # Mix with white
      r <- rgb_base[1, 1] + (1 - rgb_base[1, 1]) * fade
      g <- rgb_base[2, 1] + (1 - rgb_base[2, 1]) * fade
      b <- rgb_base[3, 1] + (1 - rgb_base[3, 1]) * fade

      color <- rgb(r, g, b)
      id <- paste(pop, pop_data$allele_group[i], sep = ".")
      result[id] <- color
    }
  }

  result
}
