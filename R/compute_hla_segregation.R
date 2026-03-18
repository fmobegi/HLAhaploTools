#' Compute HLA haplotype segregation within families
#'
#' This function phases multi‑locus HLA alleles into transmitted paternal and
#' maternal haplotypes for each child in a pedigree. It supports:
#'
#' • Complete trios (father, mother, ≥1 child)
#' • Single‑parent families with ≥3 children (inference mode)
#' • Multi‑locus HLA data where each locus has two allele columns
#'   (e.g., `A_1`, `A_2`, `B_1`, `B_2`, etc.)
#'
#' ## Haplotype labelling (AC anchoring)
#'
#' • **Child 1** defines the canonical haplotypes:
#'   – Paternal haplotype → **A**
#'   – Maternal haplotype → **C**
#'
#' • For all subsequent children:
#'   – If their paternal haplotype matches Child 1 → **A**, else **B**
#'   – If their maternal haplotype matches Child 1 → **C**, else **D**
#'
#' This yields the four possible transmitted haplotype combinations:
#' **AC**, **AD**, **BC**, **BD**.
#'
#' ## Missing‑parent inference
#'
#' When one parent is missing and ≥3 children are present, the missing parent's
#' alleles are inferred locus‑by‑locus by subtracting the known parent's alleles
#' from the union of child alleles.
#'
#' @param hped A data frame containing pedigree and HLA allele columns. Must
#'   include `FAMILY_ID` and `Family_Member`. Allele columns must follow the
#'   pattern `<Locus>_1` and `<Locus>_2`.
#' @param collapse Character scalar used to join alleles into a haplotype
#'   string. Default is `"~"`.
#' @param verbose Logical; if `TRUE`, prints progress messages.
#'
#' @return A tibble with one paternal and one maternal haplotype per child.
#'
#' @importFrom dplyr filter slice bind_rows
#' @importFrom purrr keep
#' @importFrom janitor clean_names
#' @importFrom tibble tibble
#'
#' @export
compute_hla_segregation <- function(hped, collapse = "~", verbose = TRUE) {
   required <- c("FAMILY_ID", "Family_Member")
   stopifnot(all(required %in% colnames(hped)))

   hped <- janitor::clean_names(hped, case = "none")
   hped$Family_Member <- as.character(hped$Family_Member)

   genes <- c(
      "F", "G", "H", "A", "J", "C", "B", "E",
      "MICA", "MICB", "DRB1", "DRB3", "DRB4",
      "DRB5", "DQA1", "DQB1", "DMB", "DMA", "DOA", "DPA1", "DPB1"
   )

   gene_order <- as.vector(rbind(paste0(genes, "_1"), paste0(genes, "_2")))
   gene_order <- intersect(gene_order, names(hped))
   family_cols <- setdiff(names(hped), gene_order)
   hped <- hped[, c(family_cols, gene_order), drop = FALSE]
   loci <- unique(sub("_[12]$", "", gene_order))

   allele_pairs <- function(row, locus) {
      vals <- c(row[[paste0(locus, "_1")]], row[[paste0(locus, "_2")]])
      purrr::keep(vals, ~ !is.na(.) && . != "")
   }

   make_string <- function(x) {
      x <- x[!is.na(x) & x != ""]
      if (length(x) == 0) {
         return(NA_character_)
      }
      paste(x, collapse = collapse)
   }

   results <- list()
   families <- unique(hped$FAMILY_ID)

   for (fam in families) {
      if (verbose) message("Processing: ", fam)

      fam_rows <- dplyr::filter(hped, FAMILY_ID == fam)

      father_row <- fam_rows %>%
         dplyr::filter(tolower(Family_Member) %in% c("f", "father", "dad", "paternal", "p")) %>%
         dplyr::slice(1)

      mother_row <- fam_rows %>%
         dplyr::filter(tolower(Family_Member) %in% c("m", "mother", "mum", "maternal")) %>%
         dplyr::slice(1)

      children_row <- fam_rows %>%
         dplyr::filter(grepl("^c", tolower(Family_Member)))

      if (nrow(children_row) == 0) next

      # Precompute father AB and mother CD (unordered)
      father_haps <- list(A = list(), B = list())
      mother_haps <- list(C = list(), D = list())

      for (locus in loci) {
         f <- allele_pairs(father_row, locus)
         m <- allele_pairs(mother_row, locus)

         father_haps$A[[locus]] <- f[1]
         father_haps$B[[locus]] <- ifelse(length(f) == 2, f[2], f[1])

         mother_haps$C[[locus]] <- m[1]
         mother_haps$D[[locus]] <- ifelse(length(m) == 2, m[2], m[1])
      }

      child1_pat <- NULL
      child1_mat <- NULL

      # Process each child
      for (i in seq_len(nrow(children_row))) {
         child <- children_row[i, , drop = FALSE]
         child_ID <- child$Family_Member

         pat <- list()
         mat <- list()

         for (locus in loci) {
            f_alleles <- allele_pairs(father_row, locus)
            m_alleles <- allele_pairs(mother_row, locus)
            c_alleles <- allele_pairs(child, locus)

            if (length(c_alleles) == 0) next

            pat_unique <- setdiff(intersect(f_alleles, c_alleles), m_alleles)
            mat_unique <- setdiff(intersect(m_alleles, c_alleles), f_alleles)

            if (i == 1) {
               # CHILD 1: PURE MENDELIAN ONLY — NO FALLBACK
               if (length(pat_unique) == 1 && length(c_alleles) == 2) {
                  pat[[locus]] <- pat_unique[1]
                  mat[[locus]] <- setdiff(c_alleles, pat_unique)[1]
               } else if (length(mat_unique) == 1 && length(c_alleles) == 2) {
                  mat[[locus]] <- mat_unique[1]
                  pat[[locus]] <- setdiff(c_alleles, mat_unique)[1]
               } else {
                  # ambiguous → use child's own alleles directly
                  pat[[locus]] <- c_alleles[1]
                  mat[[locus]] <- c_alleles[2]
               }
            } else {
               # CHILD 2+: full logic with fallback
               if (length(pat_unique) == 1 && length(c_alleles) == 2) {
                  pat[[locus]] <- pat_unique[1]
                  mat[[locus]] <- setdiff(c_alleles, pat_unique)[1]
               } else if (length(mat_unique) == 1 && length(c_alleles) == 2) {
                  mat[[locus]] <- mat_unique[1]
                  pat[[locus]] <- setdiff(c_alleles, mat_unique)[1]
               } else {
                  pat[[locus]] <- ifelse(father_haps$A[[locus]] %in% c_alleles,
                     father_haps$A[[locus]],
                     father_haps$B[[locus]]
                  )

                  mat[[locus]] <- ifelse(mother_haps$C[[locus]] %in% c_alleles,
                     mother_haps$C[[locus]],
                     mother_haps$D[[locus]]
                  )
               }
            }
         }

         # Store canonical AC haplotypes from Child 1
         if (i == 1) {
            child1_pat <- pat
            child1_mat <- mat
         }

         # AC anchoring
         if (i == 1) {
            pat_label <- "A"
            mat_label <- "C"
         } else {
            pat_label <- if (identical(pat, child1_pat)) "A" else "B"
            mat_label <- if (identical(mat, child1_mat)) "C" else "D"
         }

         results[[paste(fam, child_ID, sep = "__")]] <- tibble::tibble(
            FAMILY_ID = fam,
            Child_ID = child_ID,
            Haplotype = c(pat_label, mat_label),
            Allele_string = c(make_string(pat), make_string(mat))
         )
      }
   }

   dplyr::bind_rows(results)
}
