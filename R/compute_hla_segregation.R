#' Compute HLA haplotype segregation within families
#'
#' Phase HLA alleles into transmitted haplotypes for each child in a pedigree.
#' The function supports:
#' - complete trios (father, mother, children),
#' - single-parent families with three or more children (inference mode),
#' - multi-locus HLA data where each locus is represented by two allele columns
#'   (e.g., A_1, A_2).
#'
#' Haplotype labels (stable per parent across all children):
#' - **A**: paternal haplotype 1
#' - **B**: paternal haplotype 2
#' - **C**: maternal haplotype 1
#' - **D**: maternal haplotype 2
#'
#' Only the **two transmitted haplotypes** (one paternal, one maternal) are
#' returned for each child.
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
#'   - `Allele_string` (collapsed haplotype string)
#'
#' @details
#' Parental haplotypes (A/B for father, C/D for mother) are constructed once per
#' family and reused for all children. For each child, only the transmitted
#' paternal and maternal haplotypes are returned.
#'
#' @export
compute_hla_segregation <- function(hped, collapse = "~", verbose = TRUE) {
   required <- c("FAMILY_ID", "Family_Member")
   stopifnot(all(required %in% colnames(hped)))

   hped <- janitor::clean_names(hped, case = "none")

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
                     "child", "child1", "child2", "child3",
                     "c", "c1", "c2", "c3",
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

            # --- Inference mode for missing parent --------------------------------
            if (!(has_father && has_mother)) {
               if (verbose) message("\t", fam, " - Single parent detected")

               if (!enough_children) {
                  warning("\t", fam, " - Need 3+ children for inference.")
                  next
               }

               if (!has_father && has_mother) {
                  if (verbose) message("\t-> Inferring father ...")

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

                     if (length(father_alleles) >= 1) {
                        father_row[[paste0(locus, "_1")]] <- father_alleles[1]
                     }
                     if (length(father_alleles) >= 2) {
                        father_row[[paste0(locus, "_2")]] <- father_alleles[2]
                     } else if (length(father_alleles) == 1) {
                        father_row[[paste0(locus, "_2")]] <- father_alleles[1]
                     }
                  }
               } else if (has_father && !has_mother) {
                  if (verbose) message("\t-> Inferring mother ...")

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

                     if (length(mother_alleles) >= 1) {
                        mother_row[[paste0(locus, "_1")]] <- mother_alleles[1]
                     }
                     if (length(mother_alleles) >= 2) {
                        mother_row[[paste0(locus, "_2")]] <- mother_alleles[2]
                     } else if (length(mother_alleles) == 1) {
                        mother_row[[paste0(locus, "_2")]] <- mother_alleles[1]
                     }
                  }
               }
            }

            # --- PRECOMPUTE PARENTAL HAPLOTYPES (A/B for father, C/D for mother) ---
            parent_haps <- list(
               father = list(A = list(), B = list()),
               mother = list(C = list(), D = list())
            )

            for (locus in loci) {
               f_alleles <- allele_pairs(father_row, locus)
               m_alleles <- allele_pairs(mother_row, locus)

               # Father
               if (length(f_alleles) == 2) {
                  parent_haps$father$A[[locus]] <- f_alleles[1]
                  parent_haps$father$B[[locus]] <- f_alleles[2]
               } else if (length(f_alleles) == 1) {
                  parent_haps$father$A[[locus]] <- f_alleles[1]
                  parent_haps$father$B[[locus]] <- f_alleles[1]
               }

               # Mother
               if (length(m_alleles) == 2) {
                  parent_haps$mother$C[[locus]] <- m_alleles[1]
                  parent_haps$mother$D[[locus]] <- m_alleles[2]
               } else if (length(m_alleles) == 1) {
                  parent_haps$mother$C[[locus]] <- m_alleles[1]
                  parent_haps$mother$D[[locus]] <- m_alleles[1]
               }
            }

            # --- PROCESS EACH CHILD ------------------------------------------------
            for (i in seq_len(nrow(children_row))) {
               child_row <- children_row[i, , drop = FALSE]
               child_ID <- child_row$Family_Member

               child_all <- unlist(lapply(loci, function(l) allele_pairs(child_row, l)))

               patA <- unlist(parent_haps$father$A)
               patB <- unlist(parent_haps$father$B)
               matC <- unlist(parent_haps$mother$C)
               matD <- unlist(parent_haps$mother$D)

               paternal_hap <- if (all(patA %in% child_all)) {
                  "A"
               } else if (all(patB %in% child_all)) {
                  "B"
               } else {
                  NA
               }

               maternal_hap <- if (all(matC %in% child_all)) {
                  "C"
               } else if (all(matD %in% child_all)) {
                  "D"
               } else {
                  NA
               }

               # Build strings
               make_string <- function(x) {
                  if (length(x) == 0) {
                     return(NA_character_)
                  }
                  paste(x, collapse = collapse)
               }

               A_str <- make_string(patA)
               B_str <- make_string(patB)
               C_str <- make_string(matC)
               D_str <- make_string(matD)

               result_child <- tibble::tibble(
                  FAMILY_ID = fam,
                  Child_ID = child_ID,
                  Haplotype = c(paternal_hap, maternal_hap),
                  Allele_string = c(
                     if (paternal_hap == "A") A_str else B_str,
                     if (maternal_hap == "C") C_str else D_str
                  )
               ) %>% dplyr::filter(!is.na(Haplotype))

               results[[paste(fam, child_ID, sep = "__")]] <- result_child
            }

            if (verbose) message("\tFinished processing: ", fam, "\n")
         },
         error = function(e) {
            message("\t!ERROR in family ", fam, ": ", conditionMessage(e))
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
