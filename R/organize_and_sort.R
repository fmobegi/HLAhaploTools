# ~~~~~~~~~~~~~~~~~
# organize_and_sort
# ~~~~~~~~~~~~~~~~~
#' Organize and Sort HLA Data
#'
#' Internal helper function to organize columns and sort data.
#' Orders columns by gene class (Class I, Class II, non-classical, MIC)
#' and sorts records by family structure when applicable.
#'
#' @param df A data frame with cleaned and standardized HLA data
#' @param isfamilydata Logical. If TRUE, sort by family member.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A sorted and organized tibble
#' @keywords internal
#' @importFrom dplyr arrange
organize_and_sort <- function(df,
                              isfamilydata = TRUE,
                              quiet = FALSE) {
  if (!quiet) {
    message("\t📋 Organizing columns and sorting data...")
  }

  # Find ID columns based on data type
  id_cols <- if (isfamilydata) {
    # For family data, use family columns
    c("FAMILY_ID", "Family_Member")
  } else {
    # For non-family data, use any ID column we can find
    possible_id_cols <- c(
      "SampleID",
      "Sample_ID",
      "ID",
      "Id",
      "id",
      "PATIENT_ID",
      "Patient_ID",
      "Subject_ID",
      "Donor_ID",
      "donor_id",
      "RecipientID",
      "recipient_id"
    )
    # If no standard ID column, fallback to FAMILY_ID if present
    intersect(c(possible_id_cols, "FAMILY_ID"), names(df))
  }

  # Find existing allele columns
  allele_cols <- grep("_[12]$", names(df), value = TRUE)

  # Define standard gene orders within each class
  class1_order <- c("A", "B", "C")
  class2_order <- c(
    "DRB1",
    "DRB3",
    "DRB4",
    "DRB5",
    "DQA1",
    "DQB1",
    "DPA1",
    "DPB1",
    "DOA",
    "DOB",
    "DMA",
    "DMB"
  )
  nonclass_order <- c("E", "F", "G", "H", "HFE", "J", "K", "L")
  mic_order <- c("MICA", "MICB")

  # Extract gene names from columns (without _1 or _2 suffix)
  gene_names <- unique(sub("_[12]$", "", allele_cols))

  # Group genes by class
  class1_genes <- intersect(class1_order, gene_names)
  class2_genes <- intersect(class2_order, gene_names)
  nonclass_genes <- intersect(nonclass_order, gene_names)
  mic_genes <- intersect(mic_order, gene_names)
  other_genes <- setdiff(
    gene_names,
    c(class1_order, class2_order, nonclass_order, mic_order)
  )

  # Sort genes within each class according to standard order
  class1_genes <- class1_genes[match(class1_genes, class1_order)]
  class2_genes <- class2_genes[match(class2_genes, class2_order)]
  nonclass_genes <- nonclass_genes[match(nonclass_genes, nonclass_order)]
  mic_genes <- mic_genes[match(mic_genes, mic_order)]
  other_genes <- sort(other_genes) # Alphabetically sort any other genes

  # Combine all genes in order of importance
  ordered_genes <- c(
    class1_genes,
    class2_genes,
    nonclass_genes,
    mic_genes,
    other_genes
  )

  # Get columns in desired order (both _1 and _2 for each gene)
  ordered_allele_cols <- c()
  for (gene in ordered_genes) {
    col1 <- paste0(gene, "_1")
    col2 <- paste0(gene, "_2")
    if (col1 %in% names(df)) {
      ordered_allele_cols <- c(ordered_allele_cols, col1)
    }
    if (col2 %in% names(df)) {
      ordered_allele_cols <- c(ordered_allele_cols, col2)
    }
  }

  # Report gene organization if not quiet
  if (!quiet && length(gene_names) > 0) {
    message("\t   Organizing genes by HLA class:")
    if (length(class1_genes) > 0) {
      message(sprintf("\t   ▪ Class I genes: %s", paste(class1_genes, collapse = ", ")))
    }
    if (length(class2_genes) > 0) {
      message(sprintf("\t   ▪ Class II genes: %s", paste(class2_genes, collapse = ", ")))
    }
    if (length(nonclass_genes) > 0) {
      message(sprintf(
        "\t   ▪ Non-classical genes: %s",
        paste(nonclass_genes, collapse = ", ")
      ))
    }
    if (length(mic_genes) > 0) {
      message(sprintf("\t   ▪ MIC genes: %s", paste(mic_genes, collapse = ", ")))
    }
    if (length(other_genes) > 0) {
      message(sprintf("\t   ▪ Other genes: %s", paste(other_genes, collapse = ", ")))
    }
  }

  # Order columns: ID columns, then ordered allele columns, then everything else
  ordering_cols <- c(intersect(id_cols, names(df)), ordered_allele_cols)
  remaining <- setdiff(names(df), ordering_cols)

  # Apply column ordering
  df <- df[, c(ordering_cols, remaining), drop = FALSE]

  # Remove empty rows and columns
  df <- df[rowSums(is.na(df)) < ncol(df), , drop = FALSE]
  df <- df[, colSums(is.na(df)) < nrow(df), drop = FALSE]

  # Sort by family structure if applicable
  if (isfamilydata &&
    all(c("FAMILY_ID", "Family_Member") %in% names(df))) {
    # Standardize family member codes for sorting
    df$Family_Member <- sub("^[fF]$|^[fF]ather$|^[dD]ad$", "F", df$Family_Member)
    df$Family_Member <- sub("^[mM]$|^[mM]other$|^[mM]om$", "M", df$Family_Member)
    df$Family_Member <- sub("^[cC]hild[ -]?([0-9]+)$", "C\\1", df$Family_Member)
    df$Family_Member <- sub("^[cC]hild$", "C1", df$Family_Member)

    # Standard ordering: Father, Mother, Child1, Child2, ...
    child_codes <- grep("^C[0-9]+$", unique(df$Family_Member), value = TRUE)
    order_levels <- c("F", "M", sort(child_codes))

    df$Family_Member <- factor(df$Family_Member, levels = order_levels)
    df <- dplyr::arrange(df, FAMILY_ID, Family_Member)

    if (!quiet) {
      message("\t👪 Data sorted by family structure.")
    }
  }

  df
}
