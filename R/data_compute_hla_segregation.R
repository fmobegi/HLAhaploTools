#' Compute HLA haplotype segregation within families
#'
#' Performs simple rule-based HLA haplotype phasing and segregation analysis
#' across nuclear families using HPED-style genotype data. The function attempts
#' to infer parental haplotypes (A/B for father and C/D for mother) based on
#' shared alleles observed in children, then assigns the most likely haplotype
#' combination inherited by each child.
#'
#' The function supports:
#' \itemize{
#'   \item Complete trios (mother, father, children)
#'   \item Families with one missing parent (if at least 3 children are available)
#'   \item Automatic inference of missing parental alleles from offspring
#'   \item Multi-locus HLA haplotype string generation
#' }
#'
#' Expected input columns include:
#' \describe{
#'   \item{FAMILY_ID}{Unique family identifier}
#'   \item{Family_Member}{Relationship label (e.g. Mother, Father, Child1)}
#'   \item{<GENE>_1, <GENE>_2}{Diploid allele columns for each HLA locus}
#' }
#'
#' Supported loci include:
#' F, G, H, A, J, C, B, E, MICA, MICB, DRB1, DRB3, DRB4,
#' DRB5, DQA1, DQB1, DMB, DMA, DOA, DPA1, and DPB1.
#'
#' @param hped A data frame containing family structure and HLA allele calls.
#' Must contain at minimum `FAMILY_ID` and `Family_Member` columns.
#'
#' @param collapse Character separator used when concatenating alleles into
#' haplotype strings. Default is `"~"`.
#'
#' @param verbose Logical; if `TRUE`, progress messages are printed during
#' processing. Default is `TRUE`.
#'
#' @return A list containing two tibbles:
#' \describe{
#'   \item{segregation}{
#'     Child-level haplotype inheritance assignments with columns:
#'     \itemize{
#'       \item `FAMILY_ID`
#'       \item `Child_ID`
#'       \item `Haplotype`
#'       \item `Allele_string`
#'     }
#'   }
#'   \item{Allele_Strings}{
#'     Family-level phased haplotype strings for inferred paternal and maternal
#'     haplotypes.
#'   }
#' }
#'
#' @details
#' Haplotypes are assigned using the first child in the family to initialize
#' parental haplotype structure. Subsequent children are matched against the
#' four possible parental haplotype combinations:
#'
#' \itemize{
#'   \item AC
#'   \item AD
#'   \item BC
#'   \item BD
#' }
#'
#' The combination with the highest number of matching loci is selected for
#' each child.
#'
#' Families with no identified parents are skipped. Families with only one
#' parent require at least three children for missing-parent inference.
#'
#' @examples
#' \dontrun{
#' result <- compute_hla_segregation(hped)
#'
#' result$segregation
#' result$Allele_Strings
#' }
#'
#' @importFrom dplyr filter slice arrange bind_rows rowwise summarise pull
#' @importFrom purrr keep
#' @importFrom janitor clean_names
#' @importFrom stringr str_rank
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

            if (num_parents <= 1 && num_children < 3) {
               message("\t!Skipping family ", fam, ": Insufficient data for phasing.")
               next
            }

            if (nrow(children_row) == 0) next

            father_haps <- list(A = list(), B = list())
            mother_haps <- list(C = list(), D = list())

            # Infer missing parent if necessary
            if (num_parents == 1 && num_children >= 3) {
               if (verbose) {
                  message("\t-> Missing one parent, inferring from ", nrow(children_row), " children...")
               }

               known_parent_row <- if (nrow(father_row) == 1) father_row else mother_row
               missing_role <- if (nrow(father_row) == 1) "Mother" else "Father"

               inferred_row <- fam_rows[1, , drop = FALSE]
               inferred_row[1, ] <- NA
               inferred_row$FAMILY_ID <- fam
               inferred_row$Family_Member <- paste0(missing_role, " (inferred)")

               for (locus in loci) {
                  known_alleles <- allele_pairs(known_parent_row, locus)

                  all_child_alleles <- children_row %>%
                     dplyr::rowwise() %>%
                     dplyr::summarise(
                        alleles = list(allele_pairs(dplyr::cur_data(), locus)),
                        .groups = "drop"
                     ) %>%
                     dplyr::pull(alleles) %>%
                     unlist() %>%
                     unique()

                  inferred_alleles <- setdiff(all_child_alleles, known_alleles)

                  if (length(inferred_alleles) >= 1) {
                     inferred_row[[paste0(locus, "_1")]] <- inferred_alleles[1]
                  }
                  if (length(inferred_alleles) >= 2) {
                     inferred_row[[paste0(locus, "_2")]] <- inferred_alleles[2]
                  } else if (length(inferred_alleles) == 1) {
                     inferred_row[[paste0(locus, "_2")]] <- inferred_alleles[1]
                  }
               }

               if (missing_role == "Father") {
                  father_row <- inferred_row
               } else {
                  mother_row <- inferred_row
               }
            } else if (num_parents == 0) {
               warning("\t", fam, " No parents identified, cannot proceed with phasing.")
               next
            }

            # Initialize parental haplotypes using the first child
            first_child <- children_row[1, ]

            for (locus in loci) {
               child1_alleles <- sort(allele_pairs(first_child, locus))
               father_alleles <- sort(allele_pairs(father_row, locus))
               mother_alleles <- sort(allele_pairs(mother_row, locus))

               if (length(father_alleles) == 0 && length(mother_alleles) > 0) {
                  shared_m <- intersect(child1_alleles, mother_alleles)
                  mother_haps$C[[locus]] <- if (length(shared_m) > 0) shared_m[1] else mother_alleles[1]
                  remain_m <- setdiff(mother_alleles, mother_haps$C[[locus]])
                  mother_haps$D[[locus]] <- if (length(remain_m) > 0) remain_m[1] else mother_haps$C[[locus]]

                  match_idx <- match(mother_haps$C[[locus]], child1_alleles)
                  rem_child <- if (!is.na(match_idx)) child1_alleles[-match_idx] else child1_alleles
                  father_haps$A[[locus]] <- if (length(rem_child) > 0) rem_child[1] else NA_character_
                  father_haps$B[[locus]] <- NA_character_
               } else if (length(mother_alleles) == 0 && length(father_alleles) > 0) {
                  shared_f <- intersect(child1_alleles, father_alleles)
                  father_haps$A[[locus]] <- if (length(shared_f) > 0) shared_f[1] else father_alleles[1]
                  remain_f <- setdiff(father_alleles, father_haps$A[[locus]])
                  father_haps$B[[locus]] <- if (length(remain_f) > 0) remain_f[1] else father_haps$A[[locus]]

                  match_idx <- match(father_haps$A[[locus]], child1_alleles)
                  rem_child <- if (!is.na(match_idx)) child1_alleles[-match_idx] else child1_alleles
                  mother_haps$C[[locus]] <- if (length(rem_child) > 0) rem_child[1] else NA_character_
                  mother_haps$D[[locus]] <- NA_character_
               } else {
                  shared_f <- intersect(child1_alleles, father_alleles)
                  remain_f <- setdiff(father_alleles, shared_f)
                  father_haps$A[[locus]] <- if (length(shared_f) > 0) shared_f[1] else father_alleles[1]
                  father_haps$B[[locus]] <- if (length(remain_f) > 0) remain_f[1] else father_haps$A[[locus]]

                  match_idx <- match(father_haps$A[[locus]], child1_alleles)
                  child1_remain_m <- if (!is.na(match_idx)) child1_alleles[-match_idx] else child1_alleles

                  shared_m <- intersect(child1_remain_m, mother_alleles)
                  remain_m <- setdiff(mother_alleles, shared_m)
                  mother_haps$C[[locus]] <- if (length(shared_m) > 0) shared_m[1] else mother_alleles[1]
                  mother_haps$D[[locus]] <- if (length(remain_m) > 0) remain_m[1] else mother_haps$C[[locus]]
               }
            }

            fam_segregation_list <- list()

            # Score combinations for every child
            for (i in 1:nrow(children_row)) {
               child <- children_row[i, , drop = FALSE]
               child_id <- child$Family_Member

               scores <- c(AC = 0, AD = 0, BC = 0, BD = 0)
               valid_loci_count <- 0

               for (locus in loci) {
                  child_alleles <- sort(allele_pairs(child, locus))
                  if (length(child_alleles) == 0) next

                  f_A <- if (!is.null(father_haps$A[[locus]])) father_haps$A[[locus]] else NA_character_
                  f_B <- if (!is.null(father_haps$B[[locus]])) father_haps$B[[locus]] else NA_character_
                  m_C <- if (!is.null(mother_haps$C[[locus]])) mother_haps$C[[locus]] else NA_character_
                  m_D <- if (!is.null(mother_haps$D[[locus]])) mother_haps$D[[locus]] else NA_character_

                  if (is.na(f_A) && is.na(f_B)) next
                  if (is.na(m_C) && is.na(m_D)) next

                  valid_loci_count <- valid_loci_count + 1

                  combos <- list(
                     AC = sort(purrr::keep(c(f_A, m_C), ~ !is.na(.))),
                     AD = sort(purrr::keep(c(f_A, m_D), ~ !is.na(.))),
                     BC = sort(purrr::keep(c(f_B, m_C), ~ !is.na(.))),
                     BD = sort(purrr::keep(c(f_B, m_D), ~ !is.na(.)))
                  )

                  for (nm in names(combos)) {
                     if (length(child_alleles) == length(combos[[nm]]) && all(child_alleles == combos[[nm]])) {
                        scores[nm] <- scores[nm] + 1
                     }
                  }
               }

               if (sum(scores) == 0 || valid_loci_count == 0) {
                  best_combo <- "Undetermined"
                  p_label <- "Undetermined"
                  m_label <- "Undetermined"
                  p_string <- NA_character_
                  m_string <- NA_character_
               } else {
                  best_combo <- names(which.max(scores))
                  p_label <- substr(best_combo, 1, 1)
                  m_label <- substr(best_combo, 2, 2)

                  p_string <- make_string(father_haps[[p_label]])
                  m_string <- make_string(mother_haps[[m_label]])
               }

               fam_segregation_list[[child_id]] <- tibble::tibble(
                  FAMILY_ID = fam,
                  Child_ID = child_id,
                  Haplotype = c(p_label, m_label),
                  Allele_string = c(p_string, m_string)
               )
            }

            segregation_results[[fam]] <- dplyr::bind_rows(fam_segregation_list)

            allele_string_results[[fam]] <- tibble::tibble(
               FAMILY_ID = fam,
               Haplotype = c("A", "B", "C", "D"),
               Allele_String = c(
                  make_string(father_haps$A), make_string(father_haps$B),
                  make_string(mother_haps$C), make_string(mother_haps$D)
               )
            )
         },
         error = function(e) {
            message("\t!ERROR in family ", fam, ": ", conditionMessage(e))
         }
      )
   }

   return(list(
      segregation = dplyr::bind_rows(segregation_results),
      Allele_Strings = dplyr::bind_rows(allele_string_results)
   ))
}
