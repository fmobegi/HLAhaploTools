#' Compute HLA haplotype segregation within families
#'
#' This function performs two tasks: it establishes reference haplotypes (A, B, C, D)
#' based on the first child in a family, and then uses those references to phase
#' the inherited alleles for all children in the pedigree.
#'
#' The function supports:
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
#' @return A named list containing two tibbles:
#' \describe{
#'   \item{segregation}{A tibble showing the inherited haplotypes for every child:
#'     \itemize{
#'       \item \code{FAMILY_ID}: The unique family identifier.
#'       \item \code{Child_ID}: The name/ID of the child.
#'       \item \code{Haplotype}: Label of the inherited haplotype (A, B, C, or D).
#'       \item \code{Allele_string}: The specific sequence of inherited alleles.
#'     }
#'   }
#'   \item{allele_string}{A summary tibble showing the four master reference
#'     haplotypes established for each family (A, B, C, and D).}
#' }
#'
#' @details
#' The function follows a "Child 1 Reference" logic and attempts deterministic
#' phasing where possible. For ambiguous loci
#' (both parents share an allele present in the child), a simple tie-break
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

   segregation_results <- list()
   allele_string_results <- list()
   families <- unique(hped$FAMILY_ID)

   for (fam in families) {
      tryCatch(
         {
            if (verbose) message("\tProcessing: ", fam, " ...")

            fam_rows <- dplyr::filter(hped, FAMILY_ID == fam)

            father_row <- fam_rows %>%
               dplyr::filter(tolower(Family_Member) %in% c("father", "f", "dad", "d", "paternal", "p")) %>%
               dplyr::slice(1)

            mother_row <- fam_rows %>%
               dplyr::filter(tolower(Family_Member) %in% c("mother", "m", "mom", "mum", "maternal")) %>%
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

            if (nrow(children_row) == 0) next

            # Precompute father AB and mother CD (unordered)
            father_haps <- list(A = list(), B = list())
            mother_haps <- list(C = list(), D = list())

            # Infer missing parent if needed
            if (!has_father && has_mother && enough_children) {
               if (verbose) message("\t-> Missing father, inferring from ", nrow(children_row), " children")

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
                  } else if (length(father_alleles) == 1) {
                     father_row[[paste0(locus, "_2")]] <- father_alleles[1]
                  }
               }
            } else if (has_father && !has_mother && enough_children) {
               if (verbose) message("\t-> Missing mother, inferring from ", nrow(children_row), " children")

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
                  } else if (length(mother_alleles) == 1) {
                     mother_row[[paste0(locus, "_2")]] <- mother_alleles[1]
                  }
               }
            } else if (!has_father && !has_mother) {
               warning("\t", fam, " - Both parents missing, cannot proceed")
               next
            }

            if (verbose) message("\t\tFound ", nrow(children_row), " children: ", paste(children_row$Family_Member, collapse = ", "))

            for (locus in loci) {
               f <- sort(allele_pairs(father_row, locus))
               m <- sort(allele_pairs(mother_row, locus))

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
                  f_alleles <- sort(allele_pairs(father_row, locus))
                  m_alleles <- sort(allele_pairs(mother_row, locus))
                  c_alleles <- sort(allele_pairs(child, locus))

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

               segregation_results[[paste(fam, child_ID, sep = "__")]] <- tibble::tibble(
                  FAMILY_ID = fam,
                  Child_ID = child_ID,
                  Haplotype = c(pat_label, mat_label),
                  Allele_string = c(make_string(pat), make_string(mat))
               )
            }

            # Make all haplotype strings per family
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

            allele_string_results[[fam]] <- tibble::tibble(
               Allele_Strings = c(
                  A = make_string(child1_pat),
                  B = make_string(father_haps$B),
                  C = make_string(child1_mat),
                  D = make_string(mother_haps$D)
               )
            )
         },
         error = function(e) {
            message("\t!ERROR in family ", fam, ": ", conditionMessage(e))
            print(e)
         }
      )
   }
   return(list(
      segregation = dplyr::bind_rows(segregation_results),
      allele_string = dplyr::bind_rows(allele_string_results)
   ))
}
