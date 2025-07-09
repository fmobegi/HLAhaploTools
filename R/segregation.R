# ~~~~~~~~~~~~~~~~~~~~~~~
# compute_hla_segregation
# ~~~~~~~~~~~~~~~~~~~~~~~
#'
#' For each family, identify parental haplotypes, infer inheritance pattern per child,
#' and reconstruct child haplotypes directly using inferred segregation (e.g. AC, BD).
#'
#' @param hped A tibble with HLA columns and FAMILY_ID, Family_Member labels.
#' @param collapse Separator between genes in collapsed haplotype (default = "~").
#' @return A tibble with FAMILY_ID, Family_Member, Segregation, Haplotype
#'
#' @importFrom dplyr select mutate filter group_modify case_when ungroup
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer separate_rows unite pivot_wider
#'
#' @export
compute_hla_segregation <- function(hped, collapse = "~") {
  required <- c("FAMILY_ID", "Family_Member")
  stopifnot(all(required %in% colnames(hped)))

  class1 <- c("A", "B", "C")
  class2 <- c("DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1")
  nonclass <- c("F", "G", "H", "J", "E")
  mic <- c("MICA", "MICB")
  gene_order <- c(class1, class2, nonclass, mic)

  allele_cols <- grep("_[12]$", names(hped), value = TRUE)
  gene_names <- intersect(gene_order, unique(sub("_[12]$", "", allele_cols)))

  # Sort alleles within genes
  for (gene in gene_names) {
    c1 <- paste0(gene, "_1")
    c2 <- paste0(gene, "_2")
    if (all(c(c1, c2) %in% colnames(hped))) {
      sorted <- t(apply(hped[, c(c1, c2)], 1, function(r) {
        r <- r[!is.na(r) & r != "" & r != "NULL"]
        if (length(r) < 2) {
          return(c(r, NA)[1:2])
        }
        sort(r)
      }))
      hped[[c1]] <- sorted[, 1]
      hped[[c2]] <- sorted[, 2]
    }
  }

  trim_resolution <- function(allele, max_res) {
    if (is.na(allele)) {
      return(NA)
    }
    parts <- unlist(strsplit(allele, ":"))
    paste(head(parts, max_res), collapse = ":")
  }

  resolve_tag <- function(a, f1, f2, m1, m2) {
    if (is.null(a) || length(a) == 0 || is.na(a)) {
      return("x")
    }
    f_set <- na.omit(c(f1, f2))
    m_set <- na.omit(c(m1, m2))

    # Full resolution match
    if (a %in% f_set) {
      return(ifelse(a == f1, "A", ifelse(a == f2, "B", "A/B")))
    }
    if (a %in% m_set) {
      return(ifelse(a == m1, "C", ifelse(a == m2, "D", "C/D")))
    }

    max_res <- min(
      length(strsplit(a, ":")[[1]]),
      length(strsplit(f1, ":")[[1]]),
      length(strsplit(f2, ":")[[1]]),
      length(strsplit(m1, ":")[[1]]),
      length(strsplit(m2, ":")[[1]])
    )

    a_trim <- trim_resolution(a, max_res)
    f_trim <- na.omit(c(trim_resolution(f1, max_res), trim_resolution(f2, max_res)))
    m_trim <- na.omit(c(trim_resolution(m1, max_res), trim_resolution(m2, max_res)))

    # Trimmed match
    in_f <- a_trim %in% f_trim
    in_m <- a_trim %in% m_trim

    if (in_f && !in_m) {
      return("A/B")
    }
    if (!in_f && in_m) {
      return("C/D")
    }
    if (in_f && in_m) {
      return("A/B,C/D")
    }
    return(NA)
  }

  rows_out <- list()
  for (i in seq_len(nrow(hped))) {
    row <- hped[i, ]
    fam <- row$FAMILY_ID
    who <- row$Family_Member
    is_parent <- who %in% c("F", "M")

    dad <- hped[hped$FAMILY_ID == fam & hped$Family_Member == "F", ]
    mom <- hped[hped$FAMILY_ID == fam & hped$Family_Member == "M", ]

    out_row <- row
    out_row$Segregation <- ifelse(who == "F", "AB",
      ifelse(who == "M", "CD", "")
    )

    tag_row <- out_row
    tag_row$Family_Member <- "detail"

    if (!is_parent) {
      per_locus_tags <- c()
      for (gene in gene_names) {
        c1 <- row[[paste0(gene, "_1")]]
        c2 <- row[[paste0(gene, "_2")]]

        f1 <- if (paste0(gene, "_1") %in% names(dad)) dad[[paste0(gene, "_1")]] else NA
        f2 <- if (paste0(gene, "_2") %in% names(dad)) dad[[paste0(gene, "_2")]] else NA
        m1 <- if (paste0(gene, "_1") %in% names(mom)) mom[[paste0(gene, "_1")]] else NA
        m2 <- if (paste0(gene, "_2") %in% names(mom)) mom[[paste0(gene, "_2")]] else NA

        t1 <- resolve_tag(c1, f1, f2, m1, m2)
        t2 <- resolve_tag(c2, f1, f2, m1, m2)

        tag_row[[paste0(gene, "_1")]] <- t1
        tag_row[[paste0(gene, "_2")]] <- t2

        # For consensus pattern
        clean_pair <- paste0(gsub("[^ABCD]", "", t1), gsub("[^ABCD]", "", t2))
        if (nchar(clean_pair) == 2) per_locus_tags <- c(per_locus_tags, clean_pair)
      }

      sorted_patterns <- names(sort(table(per_locus_tags), decreasing = TRUE))
      valid_pattern <- NA
      for (p in sorted_patterns) {
        p_clean <- gsub("[^ABCD]", "", p)
        if (nchar(p_clean) == 2) {
          f <- substr(p_clean, 1, 1)
          m <- substr(p_clean, 2, 2)
          if (f %in% c("A", "B") && m %in% c("C", "D")) {
            valid_pattern <- p_clean
            break
          }
          if (f %in% c("C", "D") && m %in% c("A", "B")) {
            valid_pattern <- paste0(m, f)
            break
          }
        }
      }
      out_row$Segregation <- ifelse(is.na(valid_pattern), "", valid_pattern)
    }

    rows_out <- append(rows_out, list(out_row))
    if (!is_parent) rows_out <- append(rows_out, list(tag_row))
  }

  final_df <- dplyr::bind_rows(rows_out)

  df_summary_segre <- final_df |>
    dplyr::filter(Family_Member != "detail") |>
    dplyr::select(FAMILY_ID, Family_Member, Segregation)

  df_detailed_segre <- final_df

  message(sprintf(
    "\t✅ Segregation completed for %d individuals across %d families\n",
    nrow(df_summary_segre),
    length(unique(df_summary_segre$FAMILY_ID))
  ))

  df_summary_segre <- tibble::as_tibble(df_summary_segre)
  df_detailed_segre <- tibble::as_tibble(df_detailed_segre) %>%
    dplyr::select(FAMILY_ID, Family_Member, Segregation, everything())

  list(
    df_summary_segre = df_summary_segre,
    df_detailed_segre = df_detailed_segre
  )
}
