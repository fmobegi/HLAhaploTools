compute_hla_segregation <- function(hped, collapse = "~", verbose = TRUE) {
  required <- c("FAMILY_ID", "Family_Member")
  stopifnot(all(required %in% colnames(hped)))

  hped <- hped %>% janitor::clean_names(case = "none")

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

  # hped <- hped %>%
  #   dplyr::mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

  allele_pairs <- function(row, locus) {
    c(row[[paste0(locus, "_1")]], row[[paste0(locus, "_2")]]) %>%
      purrr::keep(~ !is.na(.) && . != "")
  }

  results <- list()
  families <- unique(hped$FAMILY_ID)

  for (fam in families) {
    tryCatch(
      {
        cat("\n========== PROCESSING FAMILY:", fam, "==========\n")
        fam_rows <- hped %>% filter(FAMILY_ID == fam)

        father_row <- fam_rows %>%
          filter(tolower(Family_Member) %in%
            c("father", "f", "dad", "d", "paternal", "p")) %>%
          slice(1)

        mother_row <- fam_rows %>%
          filter(tolower(Family_Member) %in%
            c("mother", "m", "mom", "mum", "maternal")) %>%
          slice(1)

        children_row <- fam_rows %>%
          filter(tolower(Family_Member) %in%
            c(
              "child", "child1", "child2", "child3", "c", "c1", "c2", "c3",
              "children", "children1", "children2", "children3"
            ))

        # Check what we have
        has_father <- nrow(father_row) > 0
        has_mother <- nrow(mother_row) > 0
        has_children <- nrow(children_row) > 0
        enough_children <- nrow(children_row) >= 3

        # Skip if no children
        if (!has_children) {
          warning(paste(fam, "- No children found"))
          next
        }

        if (has_father && has_mother) {
          if (verbose) {
            cat(fam, "- Both parents present (normal phasing)\n")
          }
        } else {
          # CASE 2: SINGLE PARENT (ONE MISSING)
          if (verbose) cat(fam, "- Single parent detected\n")

          # Check if we have enough children for inference
          if (enough_children) {
            # NESTED: Which parent is missing?
            if (!has_father && has_mother) {
              # Missing Father - Infer from Mother + Siblings
              if (verbose) {
                cat(
                  "    -> Missing father, inferring from",
                  nrow(children_row), "siblings...\n"
                )
              }

              # Create empty father row
              father_row <- fam_rows[1, ]
              father_row[1, ] <- ""
              father_row$FAMILY_ID <- fam
              father_row$Family_Member <- "Father (inferred)"

              # Inference loop for each locus
              for (locus in loci) {
                mother_alleles <- allele_pairs(mother_row, locus)

                # Collect all child alleles
                all_child_alleles <- children_row %>%
                  rowwise() %>%
                  summarise(alleles = list(allele_pairs(cur_data(), locus))) %>%
                  pull(alleles) %>%
                  unlist() %>%
                  unique()

                # Father's alleles = alleles NOT in mother
                father_alleles <- setdiff(all_child_alleles, mother_alleles)

                # Assign to father_row
                if (length(father_alleles) >= 1) {
                  father_row[[paste0(locus, "_1")]] <- father_alleles[1]
                }
                if (length(father_alleles) >= 2) {
                  father_row[[paste0(locus, "_2")]] <- father_alleles[2]
                } else if (length(father_alleles) == 1) {
                  father_row[[paste0(locus, "_2")]] <- father_alleles[1]
                }
              }

              if (verbose) cat("    -> Father inferred (Haplotype A assumed)\n")
            } else if (has_father && !has_mother) {
              # Missing Mother - Infer from Father + Siblings
              if (verbose) {
                cat(
                  "    -> Missing mother, inferring from",
                  nrow(children_row), "siblings...\n"
                )
              }

              # Create empty mother row
              mother_row <- fam_rows[1, ]
              mother_row[1, ] <- ""
              mother_row$FAMILY_ID <- fam
              mother_row$Family_Member <- "Mother (inferred)"

              # Inference loop for each locus
              for (locus in loci) {
                father_alleles <- allele_pairs(father_row, locus)

                # Collect all child alleles
                all_child_alleles <- children_row %>%
                  rowwise() %>%
                  summarise(alleles = list(allele_pairs(cur_data(), locus))) %>%
                  pull(alleles) %>%
                  unlist() %>%
                  unique()

                # Mother's alleles = alleles NOT in father
                mother_alleles <- setdiff(all_child_alleles, father_alleles)

                # Assign to mother_row
                if (length(mother_alleles) >= 1) {
                  mother_row[[paste0(locus, "_1")]] <- mother_alleles[1]
                }
                if (length(mother_alleles) >= 2) {
                  mother_row[[paste0(locus, "_2")]] <- mother_alleles[2]
                } else if (length(mother_alleles) == 1) {
                  mother_row[[paste0(locus, "_2")]] <- mother_alleles[1]
                }
              }

              if (verbose) cat("    -> Mother inferred (Haplotype C assumed)\n")
            } else {
              # Both parents missing
              warning(paste(fam, "- Both parents missing, cannot proceed"))
              next
            }
          } else {
            # Single parent but NOT enough children
            warning(paste(
              fam, "- Single parent with only", nrow(children_row),
              "child(ren). Need 3+ for inference."
            ))
            next
          }
        }

        num_children <- nrow(children_row)
        cat(
          " Found", num_children, "children:",
          paste(children_row$Family_Member, collapse = ", "), "\n"
        )

        for (i in seq_len(nrow(children_row))) {
          child_row <- children_row[i, ]
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

            if (length(pat_unique) > 0) {
              pat_allele <- pat_unique[1]
            }

            if (length(mat_unique) > 0) {
              mat_allele <- mat_unique[1]
            }

            if (!is.na(pat_allele) && is.na(mat_allele)) {
              remaining <- setdiff(c_alleles, pat_allele)
              if (length(remaining) > 0) mat_allele <- remaining[1]
            }

            if (is.na(pat_allele) && !is.na(mat_allele)) {
              remaining <- setdiff(c_alleles, mat_allele)
              if (length(remaining) > 0) pat_allele <- remaining[1]
            }

            # Handle ambiguous cases
            if (is.na(pat_allele) && is.na(mat_allele)) {
              if (length(ambiguous) >= 1) {
                pat_allele <- ambiguous[1]
              }
              if (length(ambiguous) >= 2) {
                mat_allele <- ambiguous[2]
              } else if (length(c_alleles) >= 2) {
                mat_allele <- c_alleles[2]
              } else {
                mat_allele <- c_alleles[1]
              }
            }

            # Store in Haplotype list
            if (!is.na(pat_allele) && pat_allele != "") {
              hapA[[locus]] <- pat_allele
            }
            if (!is.na(mat_allele) && mat_allele != "") {
              hapC[[locus]] <- mat_allele
            }

            # Non-transmitted alleles (WITH THE FIX)
            if (length(f_alleles) > 0) {
              pat_non <- setdiff(f_alleles, pat_allele)
              if (length(pat_non) > 0) {
                hapB[[locus]] <- pat_non[1]
              } else {
                # Homozygous case: assign the same allele
                hapB[[locus]] <- pat_allele
              }
            }

            if (length(m_alleles) > 0) {
              mat_non <- setdiff(m_alleles, mat_allele)
              if (length(mat_non) > 0) {
                hapD[[locus]] <- mat_non[1]
              } else {
                # Homozygous case: assign the same allele
                hapD[[locus]] <- mat_allele
              }
            }
          }

          # Create haplotype strings
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

          # Create haplotype strings in chromosomal order
          A_str <- make_string(hapA)
          B_str <- make_string(hapB)
          C_str <- make_string(hapC)
          D_str <- make_string(hapD)

          # Create output dataframe for the child
          result_child <- tibble::tibble(
            FAMILY_ID = fam,
            Child_ID = child_ID,
            Haplotype = c("A", "B", "C", "D"),
            Allele_string = c(A_str, B_str, C_str, D_str)
          )

          # Store results
          results[[paste(fam, child_ID, sep = "__")]] <- result_child
        }

        if (verbose) cat("\n")
        cat("========== FINISHED FAMILY:", fam, "==========\n")
      },
      error = function(e) {
        cat("!!!!! ERROR in family", fam, ":", conditionMessage(e), "\n")
        print(e)
      }
    )
  }

  final <- dplyr::bind_rows(results)

  if (nrow(final) == 0) {
    message("No haplotypes produced")
    return(NULL)
  }

  return(final)
}
