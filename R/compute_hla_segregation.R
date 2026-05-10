#' Compute HLA haplotype segregation within families
#'
#' This function performs multi-locus HLA phasing within pedigrees. It establishes
#' master reference haplotypes (A, B, C, D) based on the first child in each
#' family and then uses a linkage-aware scoring system to phase inherited
#' alleles for all siblings.
#'
#' @description
#' The function supports:
#' \itemize{
#'   \item Complete trios (Father, Mother, Children).
#'   \item Single-parent families with three or more children (Inference Mode).
#'   \item Multi-locus data where each locus has two allele columns (e.g., A_1, A_2).
#' }
#'
#' Haplotype definitions:
#' \itemize{
#'   \item \strong{A}: Paternal haplotype transmitted to Child 1.
#'   \item \strong{B}: Paternal haplotype NOT transmitted to Child 1.
#'   \item \strong{C}: Maternal haplotype transmitted to Child 1.
#'   \item \strong{D}: Maternal haplotype NOT transmitted to Child 1.
#' }
#'
#' @param hped A data frame containing pedigree information and HLA alleles.
#'   Must include \code{FAMILY_ID} and \code{Family_Member}. Allele columns
#'   must follow the pattern \code{<Locus>_1} and \code{<Locus>_2}.
#' @param collapse Character scalar used to join alleles into a haplotype
#'   string (default: \code{"~"}).
#' @param verbose Logical; if \code{TRUE}, prints progress messages to the
#'   console (default: \code{TRUE}).
#'
#' @return A named list containing two tibbles:
#' \describe{
#'   \item{segregation}{A tibble with phased results for every child:
#'     \itemize{
#'       \item \code{FAMILY_ID}: The unique family identifier.
#'       \item \code{Child_ID}: The ID/Member name of the child.
#'       \item \code{Haplotype}: The assigned haplotype label (A, B, C, or D).
#'       \item \code{Allele_strings}: The collapsed sequence of alleles for that haplotype.
#'     }
#'   }
#'   \item{allele_string}{A dictionary tibble showing the master reference
#'     haplotypes (A, B, C, D) established for each family.}
#' }
#'
#' @details
#' The function utilizes a "Child 1 Reference" logic. For families missing a
#' parent, the missing alleles are inferred by subtracting the known parent's
#' alleles from the collective pool of all children's alleles. Phasing for
#' siblings is determined by calculating the best-fit match score across all
#' four possible Mendelian combinations (AC, AD, BC, BD), which provides
#' robustness against missing data or homozygous loci.
#'
#' @examples
#' \dontrun{
#' # hped must contain FAMILY_ID, Family_Member, and locus_1/locus_2 columns
#' res <- compute_hla_segregation(hped, collapse = "~", verbose = TRUE)
#'
#' # Access results
#' phased_children <- res$segregation
#' family_haplotypes <- res$allele_string
#' }
#'
#' @importFrom dplyr filter slice rowwise summarise pull bind_rows rename_with matches select contains
#' @importFrom purrr keep
#' @importFrom janitor clean_names
#' @importFrom tibble tibble
#' @importFrom stringr str_rank
#' @export
compute_hla_segregation <- function(hped, collapse = "~", verbose = TRUE) {
   required <- dplyr::select(
      hped, dplyr::contains("^family_id$"), dplyr::contains("^family_member$")
   )
   stopifnot(all(required %in% colnames(hped)))

   hped <- janitor::clean_names(hped, case = "none")
   hped$Family_Member <- as.character(hped$Family_Member)

   hped <- hped %>%
      dplyr::rename_with(~"FAMILY_ID", dplyr::matches("^family_id$")) %>%
      dplyr::rename_with(~"Family_Member", dplyr::matches("^family_member$"))

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
      prefixed_alleles <- sapply(loci_present, function(loc) {
         val <- hap_list[[loc]]
         if (is.na(val) || val == "") {
            return(NULL)
         }
         paste0(loc, "*", val)
      })

      prefixed_alleles <- unlist(prefixed_alleles)

      if (length(prefixed_alleles) == 0) {
         return(NA_character_)
      }

      paste(prefixed_alleles, collapse = collapse)
   }

   segregation_results <- list()
   allele_string_results <- list()
   families <- unique(hped$FAMILY_ID)

   for (fam in families) {
      tryCatch(
         {
            if (verbose) message("\tProcessing: ", fam, " ...")

            fam_rows <- dplyr::filter(hped, FAMILY_ID == fam)

            mother_row <- fam_rows %>%
               dplyr::filter(grepl(
                  "^(mother|m|mom|mum|maternal)$",
                  tolower(Family_Member)
               )) %>%
               dplyr::slice(1)

            father_row <- fam_rows %>%
               dplyr::filter(grepl(
                  "^(father|f|dad|paternal|p)$",
                  tolower(Family_Member)
               )) %>%
               dplyr::slice(1)

            children_row <- fam_rows %>%
               dplyr::filter(grepl(
                  "^(child(ren)?|c)[0-9]*$",
                  tolower(Family_Member)
               )) %>%
               dplyr::arrange(stringr::str_rank(Family_Member, numeric = TRUE))

            num_parents <- nrow(mother_row) + nrow(father_row)
            num_children <- nrow(children_row)

            # <1 parent and <3 children
            if (num_parents <= 1 && num_children < 3) {
               message("\t!Skipping family ", fam, ": Insufficient data for phasing.")
               next
            }

            if (nrow(children_row) == 0) next # no children

            # Precompute father AB and mother CD (unordered)
            father_haps <- list(A = list(), B = list())
            mother_haps <- list(C = list(), D = list())

            # Inferring missing parent
            if (num_parents == 1 && num_children >= 3) {
               if (verbose) {
                  message(
                     "\t-> Missing one parent, inferring from ",
                     nrow(children_row), " children..."
                  )
               }

               known_parent_row <- if (nrow(father_row) == 1) father_row else mother_row
               missing_role <- if (nrow(father_row) == 1) "Mother" else "Father"

               # Create the structure
               inferred_row <- fam_rows[1, , drop = FALSE]
               inferred_row[1, ] <- NA
               inferred_row$FAMILY_ID <- fam
               inferred_row$Family_Member <- paste0(missing_role, " (inferred)")

               for (locus in loci) {
                  known_alleles <- allele_pairs(known_parent_row, locus)

                  # Get every allele present in any child at this locus
                  all_child_alleles <- children_row %>%
                     dplyr::rowwise() %>%
                     dplyr::summarise(
                        alleles = list(allele_pairs(cur_data(), locus)),
                        .groups = "drop"
                     ) %>%
                     dplyr::pull(alleles) %>%
                     unlist() %>%
                     unique()

                  # Missing Parent = (All Child Alleles) minus (Known Parent Alleles)
                  inferred_alleles <- setdiff(all_child_alleles, known_alleles)

                  if (length(
                     inferred_alleles
                  ) >= 1) {
                     inferred_row[[paste0(locus, "_1")]] <-
                        inferred_alleles[1]
                  }

                  if (length(inferred_alleles) >= 2) {
                     inferred_row[[paste0(locus, "_2")]] <- inferred_alleles[2]
                  } else if (length(inferred_alleles) == 1) {
                     inferred_row[[paste0(locus, "_2")]] <- inferred_alleles[1]
                  }
               }

               if (missing_role == "Father") {
                  father_row <-
                     inferred_row
               } else {
                  mother_row <- inferred_row
               }
            } else if (num_parents == 0) {
               warning("\t", fam, "No parents identified, cannot proceed with phasing.")
               next
            }

            first_child <- children_row[1, ]

            for (locus in loci) {
               child1_alleles <- sort(allele_pairs(first_child, locus))
               father_alleles <- sort(allele_pairs(father_row, locus))
               mother_alleles <- sort(allele_pairs(mother_row, locus))

               # CASE 1: Father missing only
               if (length(father_alleles) == 0 && length(mother_alleles) > 0) {
                  shared_m <- intersect(child1_alleles, mother_alleles)
                  mother_haps$C[[locus]] <-
                     if (length(shared_m) > 0) shared_m[1] else mother_alleles[1]
                  remain_m <- setdiff(mother_alleles, mother_haps$C[[locus]])
                  mother_haps$D[[locus]] <-
                     if (length(remain_m) > 0) remain_m[1] else mother_haps$C[[locus]]

                  match_idx <- match(mother_haps$C[[locus]], child1_alleles)
                  rem_child <-
                     if (!is.na(match_idx)) child1_alleles[-match_idx] else child1_alleles
                  father_haps$A[[locus]] <-
                     if (length(rem_child) > 0) rem_child[1] else NA_character_
                  father_haps$B[[locus]] <- NA_character_

                  # CASE 2: Mother missing only
               } else if (length(mother_alleles) == 0 && length(father_alleles) > 0) {
                  shared_f <- intersect(child1_alleles, father_alleles)
                  father_haps$A[[locus]] <-
                     if (length(shared_f) > 0) shared_f[1] else father_alleles[1]
                  remain_f <- setdiff(father_alleles, father_haps$A[[locus]])
                  father_haps$B[[locus]] <-
                     if (length(remain_f) > 0) remain_f[1] else father_haps$A[[locus]]

                  match_idx <- match(father_haps$A[[locus]], child1_alleles)
                  rem_child <-
                     if (!is.na(match_idx)) child1_alleles[-match_idx] else child1_alleles
                  mother_haps$C[[locus]] <-
                     if (length(rem_child) > 0) rem_child[1] else NA_character_
                  mother_haps$D[[locus]] <- NA_character_

                  # CASE 3: normal
               } else {
                  shared_f <- intersect(child1_alleles, father_alleles)
                  remain_f <- setdiff(father_alleles, shared_f)
                  father_haps$A[[locus]] <- shared_f[1]
                  father_haps$B[[locus]] <- ifelse(length(remain_f) > 0,
                     remain_f[1], father_haps$A[[locus]]
                  )

                  # If mother shares same allele as father & child
                  match_idx <- match(father_haps$A[[locus]], child1_alleles)
                  child1_remain_m <-
                     if (!is.na(match_idx)) child1_alleles[-match_idx] else child1_alleles

                  shared_m <- intersect(child1_remain_m, mother_alleles)
                  remain_m <- setdiff(mother_alleles, shared_m)
                  mother_haps$C[[locus]] <- shared_m[1]
                  mother_haps$D[[locus]] <-
                     ifelse(length(remain_m) > 0,
                        remain_m[1], mother_haps$C[[locus]]
                     )
               }
            }

            fam_segregation_list <- list()

            for (i in 1:nrow(children_row)) {
               child <- children_row[i, , drop = FALSE]
               child_id <- child$Family_Member

               scores <- c(AC = 0, AD = 0, BC = 0, BD = 0)

               for (locus in loci) {
                  child_alleles <- sort(allele_pairs(child, locus))
                  if (length(child_alleles) == 0) next

                  # Define all combinations
                  combos <- list(
                     AC = sort(c(father_haps$A[[locus]], mother_haps$C[[locus]])),
                     AD = sort(c(father_haps$A[[locus]], mother_haps$D[[locus]])),
                     BC = sort(c(father_haps$B[[locus]], mother_haps$C[[locus]])),
                     BD = sort(c(father_haps$B[[locus]], mother_haps$D[[locus]]))
                  )

                  for (nm in names(combos)) {
                     if (identical(child_alleles, combos[[nm]])) {
                        scores[nm] <- scores[nm] + 1
                     }
                  }
               }

               best_combo <- names(which.max(scores))
               p_label <- substr(best_combo, 1, 1)
               m_label <- substr(best_combo, 2, 2)

               p_string <- make_string(father_haps[[p_label]])
               m_string <- make_string(mother_haps[[m_label]])

               fam_segregation_list[[child_id]] <- tibble::tibble(
                  FAMILY_ID = fam,
                  Child_ID = child_id,
                  Haplotype = c(p_label, m_label),
                  Allele_Strings = c(p_string, m_string)
               )
            }

            segregation_results[[fam]] <- dplyr::bind_rows(fam_segregation_list)

            allele_string_results[[fam]] <- tibble::tibble(
               Allele_Strings = c(
                  make_string(father_haps$A), make_string(father_haps$B),
                  make_string(mother_haps$C), make_string(mother_haps$D)
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
      Segregation = dplyr::bind_rows(segregation_results),
      Allele_Strings = dplyr::bind_rows(allele_string_results)
   ))
}
