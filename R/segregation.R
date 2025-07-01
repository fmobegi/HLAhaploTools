# ~~~~~~~~~~~~~~~~~~~~~~~
# compute_hla_segregation
# ~~~~~~~~~~~~~~~~~~~~~~~
#' Compute Sorted and Cleaned HLA Haplotypes (Verbose)
#'
#' Extracts haplotypes by sorting allele pairs, decoding multi-alleles,
#' excluding empty/nulls, and collapsing by gene group within families.
#'
#' @param hped A data frame with paired HLA columns, including FAMILY_ID and Family_Member.
#' @param collapse Separator between genes in collapsed haplotype (default = "~").
#'
#' @return A tibble with FAMILY_ID, Family_Member, haplotype_id, haplotype, and per-locus values.
#' @export
#'
#' @importFrom dplyr select mutate filter all_of bind_rows
#' @importFrom tidyr pivot_longer separate_rows unite pivot_wider
compute_hla_segregation <- function(hped, collapse = "~") {
  message("\t🔄 START Family segregation Analysis")

  required <- c("FAMILY_ID", "Family_Member")
  stopifnot(all(required %in% colnames(hped)))

  class1 <- c("A", "B", "C")
  class2 <- c("DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1")
  nonclass <- c("F", "G", "H", "J", "E")
  mic <- c("MICA", "MICB")
  gene_order <- c(class1, class2, nonclass, mic)

  allele_cols <- grep("_[12]$", names(hped), value = TRUE)
  gene_names <- intersect(gene_order, unique(sub("_[12]$", "", allele_cols)))

  message("\t🔬 1. Sorting allele_1 and allele_2 within each gene...")
  for (gene in gene_names) {
    col1 <- paste0(gene, "_1")
    col2 <- paste0(gene, "_2")
    if (all(c(col1, col2) %in% names(hped))) {
      sorted <- t(apply(hped[, c(col1, col2)], 1, function(r) {
        r <- r[!is.na(r) & r != "" & r != "NULL"]
        if (length(r) < 2) {
          return(c(NA, NA))
        }
        sort(r)
      }))
      hped[[col1]] <- sorted[, 1]
      hped[[col2]] <- sorted[, 2]
    }
  }

  message("\t📐 2. Reshaping to long format and decoding MAC strings...")
  long_df <- hped %>%
    dplyr::select(FAMILY_ID, Family_Member, dplyr::all_of(allele_cols)) %>%
    tidyr::pivot_longer(
      -c(FAMILY_ID, Family_Member),
      names_to   = c("locus", "copy"),
      names_sep  = "_",
      values_to  = "allele"
    ) %>%
    tidyr::separate_rows(allele, sep = "/") %>%
    dplyr::filter(!is.na(allele) &
      allele != "" & allele != "NULL") %>%
    dplyr::mutate(copy = paste0("H", copy)) %>%
    tidyr::unite("locus_copy", locus, copy, sep = "_") %>%
    tidyr::pivot_wider(names_from = locus_copy, values_from = allele)

  message("\t🧬 3. Constructing phased haplotypes by collapsing loci...")
  get_cols <- function(copy) {
    cols <- paste0(gene_names, "_", copy)
    intersect(cols, names(long_df))
  }

  hap1_cols <- get_cols("H1")
  hap2_cols <- get_cols("H2")

  clean_and_unite <- function(df, cols, label) {
    df %>%
      dplyr::select(FAMILY_ID, Family_Member, dplyr::all_of(cols)) %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(cols), ~ ifelse(
        . %in% c("NA", "NULL", "", NA), NA_character_, .
      ))) %>%
      tidyr::unite("haplotype",
        dplyr::all_of(cols),
        sep = collapse,
        na.rm = TRUE
      ) %>%
      dplyr::mutate(
        haplotype_id = paste(
          FAMILY_ID,
          Family_Member,
          label,
          sep = "_"
        )
      )
  }

  hap1 <- clean_and_unite(long_df, hap1_cols, "H1")
  hap2 <- clean_and_unite(long_df, hap2_cols, "H2")

  message("\t✅ END Family segregation: Merging results...")
  result <- dplyr::bind_rows(hap1, hap2) %>%
    dplyr::select(
      FAMILY_ID,
      Family_Member,
      haplotype_id,
      haplotype,
      dplyr::everything()
    )

  message(sprintf(
    "\t🏁 Extracted %d haplotypes across %d records.\n",
    nrow(result),
    nrow(hped)
  ))
  return(result)
}
