#' Compute HLA haplotype segregation within families
#'
#' Phase HLA alleles into transmitted and non-transmitted haplotypes for each
#' child in a pedigree. The function supports:
#' - complete trios (father, mother, children),
#' - single-parent families with three or more children (inference mode),
#' - multi-locus HLA data where each locus is represented by two allele columns
#'   (e.g., A_1, A_2).
#'
#' Haplotype labels:
#' - **A**: paternal transmitted
#' - **B**: paternal non-transmitted
#' - **C**: maternal transmitted
#' - **D**: maternal non-transmitted
#'
#' When one parent is missing and there are at least three children, the missing
#' parent's alleles are inferred by subtracting the known parent's alleles from
#' the union of child alleles at each locus.
#'
#' @param hped A data frame containing pedigree and HLA allele columns. Must
#'   include `FAMILY_ID` and `Family_Member`. Allele columns should follow the
#'   pattern `<Locus>_1` and `<Locus>_2`.
#' @param collapse Character scalar used to join alleles into a haplotype
#'   string (default: `"~"`).
#' @param verbose Logical; if `TRUE` prints progress messages (default: `TRUE`).
#'
#' @return A tibble with columns:
#'   - `FAMILY_ID`
#'   - `Child_ID`
#'   - `Haplotype` (one of `A`, `B`, `C`, `D`)
#'   - `Allele_string` (collapsed haplotype string or `NA` if no alleles)
#'
#' @details
#' The function attempts deterministic phasing where possible. For ambiguous
#' loci (both parents share an allele present in the child), a simple tie-break
#' heuristic is used to assign transmitted alleles. Non-transmitted alleles are
#' taken from the parent's remaining allele; in homozygous cases the same allele
#' may be used for both transmitted and non-transmitted haplotypes.
#'
#' @examples
#' \dontrun{
#' # hped must contain FAMILY_ID, Family_Member and locus_1 / locus_2 columns
#' res <- compute_hla_segregation(hped, collapse = "-", verbose = TRUE)
#' }
#'
#' @importFrom cli cli_alert_info cli_alert_success cli_warn
#' @importFrom dplyr filter slice rowwise summarise pull bind_rows mutate across everything
#' @importFrom purrr keep
#' @importFrom janitor clean_names
#' @importFrom tibble tibble
#' @export
compute_hla_segregation <- function(hped, collapse = "~", verbose = TRUE) {
   required <- c("FAMILY_ID", "Family_Member")
   stopifnot(all(required %in% colnames(hped)))

   # Keep original column names but clean them for safe indexing
   hped <- janitor::clean_names(hped, case = "none")

   # Preferred gene order (keeps compatibility with typical HLA panels)
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

   # Helper to extract allele pair for a locus from a row (vector or tibble row)
   allele_pairs <- function(row, locus) {
      vals <- c(row[[paste0(locus, "_1")]], row[[paste0(locus, "_2")]])
      purrr::keep(vals, ~ !is.na(.) && . != "")
   }

   results <- list()
   families <- unique(hped$FAMILY_ID)

   for (fam in families) {
      tryCatch(
         {
            if (verbose) message("\tProcessing: ", fam, " ...")

            fam_rows <- dplyr::filter(hped, FAMILY_ID == fam)

            father_row <- fam_rows %>%
               dplyr::filter(tolower(Family_Member) %in%
                  c("father", "f", "dad", "d", "paternal", "p")) %>%
               dplyr::slice(1)

            mother_row <- fam_rows %>%
               dplyr::filter(tolower(Family_Member) %in%
                  c("mother", "m", "mom", "mum", "maternal")) %>%
               dplyr::slice(1)

            children_row <- fam_rows %>%
               dplyr::filter(tolower(Family_Member) %in%
                  c(
                     "child", "child1", "child2", "child3", "c", "c1", "c2", "c3",
                     "children", "children1", "children2", "children3"
                  ))

            has_father <- nrow(father_row) > 0
            has_mother <- nrow(mother_row) > 0
            has_children <- nrow(children_row) > 0
            enough_children <- nrow(children_row) >= 3

            if (!has_children) {
               warning("\t", fam, " - No children found")
               next
            }

            if (has_father && has_mother) {
               if (verbose) message("\t", fam, " - Both parents present (normal phasing)")
            } else {
               if (verbose) message("\t", fam, " - Single parent detected")

               if (!enough_children) {
                  warning(
                     "\t", fam, " - Single parent with only ", nrow(children_row),
                     " child(ren). Need 3+ for inference."
                  )
                  next
               }

               # Infer missing parent from children + known parent
               if (!has_father && has_mother) {
                  if (verbose) {
                     message(
                        "\t-> Missing father, inferring from ",
                        nrow(children_row), " siblings ..."
                     )
                  }

                  father_row <- fam_rows[1, , drop = FALSE]
                  father_row[1, ] <- ""
                  father_row$FAMILY_ID <- fam
                  father_row$Family_Member <- "Father (inferred)"

                  for (locus in loci) {
                     mother_alleles <- allele_pairs(mother_row, locus)

                     all_child_alleles <- children_row %>%
                        dplyr::rowwise() %>%
                        dplyr::summarise(alleles = list(allele_pairs(cur_data(), locus))) %>%
                        dplyr::pull(alleles) %>%
                        unlist() %>%
                        unique()

                     father_alleles <- setdiff(all_child_alleles, mother_alleles)

                     if (length(father_alleles) >= 1) father_row[[paste0(locus, "_1")]] <- father_alleles[1]
                     if (length(father_alleles) >= 2) {
                        father_row[[paste0(locus, "_2")]] <- father_alleles[2]
                     } else if (length(father_alleles) == 1) father_row[[paste0(locus, "_2")]] <- father_alleles[1]
                  }

                  if (verbose) message("\t-> Father inferred (haplotype A assumed)")
               } else if (has_father && !has_mother) {
                  if (verbose) {
                     message(
                        "\t-> Missing mother, inferring from ",
                        nrow(children_row), " siblings ..."
                     )
                  }

                  mother_row <- fam_rows[1, , drop = FALSE]
                  mother_row[1, ] <- ""
                  mother_row$FAMILY_ID <- fam
                  mother_row$Family_Member <- "Mother (inferred)"

                  for (locus in loci) {
                     father_alleles <- allele_pairs(father_row, locus)

                     all_child_alleles <- children_row %>%
                        dplyr::rowwise() %>%
                        dplyr::summarise(alleles = list(allele_pairs(cur_data(), locus))) %>%
                        dplyr::pull(alleles) %>%
                        unlist() %>%
                        unique()

                     mother_alleles <- setdiff(all_child_alleles, father_alleles)

                     if (length(mother_alleles) >= 1) mother_row[[paste0(locus, "_1")]] <- mother_alleles[1]
                     if (length(mother_alleles) >= 2) {
                        mother_row[[paste0(locus, "_2")]] <- mother_alleles[2]
                     } else if (length(mother_alleles) == 1) mother_row[[paste0(locus, "_2")]] <- mother_alleles[1]
                  }

                  if (verbose) message("\t-> Mother inferred (haplotype C assumed)")
               } else {
                  warning("\t", fam, " - Both parents missing, cannot proceed")
                  next
               }
            }

            num_children <- nrow(children_row)
            if (verbose) {
               message(
                  "\t\tFound ", num_children, " children: ",
                  paste(children_row$Family_Member, collapse = ", ")
               )
            }

            for (i in seq_len(nrow(children_row))) {
               child_row <- children_row[i, , drop = FALSE]
               child_ID <- child_row$Family_Member

               hapA <- list()
               hapB <- list()
               hapC <- list()
               hapD <- list()

               for (locus in loci) {
                  f_alleles <- allele_pairs(father_row, locus)
                  m_alleles <- allele_pairs(mother_row, locus)
                  c_alleles <- allele_pairs(child_row, locus)

                  if (length(c_alleles) == 0) next

                  pat_trans <- intersect(f_alleles, c_alleles)
                  mat_trans <- intersect(m_alleles, c_alleles)

                  pat_unique <- setdiff(pat_trans, mat_trans)
                  mat_unique <- setdiff(mat_trans, pat_trans)
                  ambiguous <- intersect(pat_trans, mat_trans)

                  pat_allele <- NA_character_
                  mat_allele <- NA_character_

                  if (length(pat_unique) > 0) pat_allele <- pat_unique[1]
                  if (length(mat_unique) > 0) mat_allele <- mat_unique[1]

                  if (!is.na(pat_allele) && is.na(mat_allele)) {
                     remaining <- setdiff(c_alleles, pat_allele)
                     if (length(remaining) > 0) mat_allele <- remaining[1]
                  }

                  if (is.na(pat_allele) && !is.na(mat_allele)) {
                     remaining <- setdiff(c_alleles, mat_allele)
                     if (length(remaining) > 0) pat_allele <- remaining[1]
                  }

                  # Ambiguous assignment fallback
                  if (is.na(pat_allele) && is.na(mat_allele)) {
                     if (length(ambiguous) >= 1) pat_allele <- ambiguous[1]
                     if (length(ambiguous) >= 2) {
                        mat_allele <- ambiguous[2]
                     } else if (length(c_alleles) >= 2) {
                        mat_allele <- c_alleles[2]
                     } else {
                        mat_allele <- c_alleles[1]
                     }
                  }

                  if (!is.na(pat_allele) && pat_allele != "") hapA[[locus]] <- pat_allele
                  if (!is.na(mat_allele) && mat_allele != "") hapC[[locus]] <- mat_allele

                  # Non-transmitted alleles
                  if (length(f_alleles) > 0) {
                     pat_non <- setdiff(f_alleles, pat_allele)
                     if (length(pat_non) > 0) {
                        hapB[[locus]] <- pat_non[1]
                     } else {
                        hapB[[locus]] <- pat_allele
                     }
                  }

                  if (length(m_alleles) > 0) {
                     mat_non <- setdiff(m_alleles, mat_allele)
                     if (length(mat_non) > 0) {
                        hapD[[locus]] <- mat_non[1]
                     } else {
                        hapD[[locus]] <- mat_allele
                     }
                  }
               }

               make_string <- function(hap_list) {
                  if (length(hap_list) == 0) {
                     return(NA_character_)
                  }
                  loci_present <- intersect(loci, names(hap_list))
                  alleles <- unname(unlist(hap_list[loci_present]))
                  alleles <- alleles[!is.na(alleles) & alleles != ""]
                  if (length(alleles) == 0) {
                     return(NA_character_)
                  }
                  paste(alleles, collapse = collapse)
               }

               A_str <- make_string(hapA)
               B_str <- make_string(hapB)
               C_str <- make_string(hapC)
               D_str <- make_string(hapD)

               result_child <- tibble::tibble(
                  FAMILY_ID = fam,
                  Child_ID = child_ID,
                  Haplotype = c("A", "B", "C", "D"),
                  Allele_string = c(A_str, B_str, C_str, D_str)
               )

               results[[paste(fam, child_ID, sep = "__")]] <- result_child
            }

            if (verbose) message("\tFinished processing: ", fam, "\n")
         },
         error = function(e) {
            message("\t!ERROR in family ", fam, ": ", conditionMessage(e))
            print(e)
         }
      )
   }

   final <- dplyr::bind_rows(results)

   if (nrow(final) == 0) {
      warning("No haplotypes produced")
      return(NULL)
   }

   return(final)
}
