# ~~~~~~~~~~~~~~~~~~~~
# fill_missing_alleles
# ~~~~~~~~~~~~~~~~~~~~
#' Fill Missing Alleles with DRB1 vs DRB3/4/5 Compliance Checks
#'
#' Internal helper function to:
#' - Copy first allele to second allele position when missing,
#' - Reorder DRB alleles to ensure DRBx_1 is populated when present only in _2,
#' - Fill DRB3/4/5 alleles based on DRB1 linkage rules,
#' - Identify DRB1-DRB3/4/5 inconsistencies.
#'
#' @param df A data frame with standardized allele columns (e.g., A_1, DRB1_2, etc.).
#' @param id_cols Optional character vector of column names to identify samples (e.g., "SampleID").
#'                If NULL, common ID column names will be automatically detected.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{cleaned}{The updated tibble with alleles filled and reordered.}
#'   \item{qc}{A tibble with ID columns and a `DRB145_Compliant` logical flag.}
#' }
#' @keywords internal
fill_missing_alleles <- function(df, id_cols = NULL, quiet = FALSE) {
  if (!quiet) message("\t🔄 Starting allele fill with DRB reordering and compliance check...")

  # Auto-detect ID columns if not provided
  if (is.null(id_cols)) {
    id_candidates <- c(
      "FAMILY_ID", "Family_ID", "family_id", "FamilyID", "familyid", "Family", "family",
      "Family_Member", "FamilyMember", "Member", "MEMBER", "Relationship", "RELATIONSHIP",
      "SampleID", "Sample_ID", "ID", "Id", "id", "PATIENT_ID", "Patient_ID", "Subject_ID",
      "Donor_ID", "donor_id", "RecipientID", "recipient_id"
    )
    id_cols <- intersect(id_candidates, names(df))
  }

  # Helper function to get sample identifier string
  get_sample_id <- function(i) {
    if (is.null(id_cols)) {
      return(paste0("Row", i))
    }
    id_vals <- sapply(id_cols, function(col) {
      val <- df[[col]][i]
      if (is.na(val) || val == "") {
        return(NULL)
      }
      as.character(val)
    })
    id_vals <- id_vals[!sapply(id_vals, is.null)]
    if (length(id_vals) == 0) {
      return(paste0("Row", i))
    }
    paste(id_vals, collapse = "_")
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Step 1: Reorder DRB alleles
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  reorder_drb_alleles <- function(df_local, quiet_local) {
    drb_genes <- c("DRB1", "DRB3", "DRB4", "DRB5")
    swapped_total <- 0

    for (gene in drb_genes) {
      a1 <- paste0(gene, "_1")
      a2 <- paste0(gene, "_2")
      if (all(c(a1, a2) %in% names(df_local))) {
        for (i in seq_len(nrow(df_local))) {
          if ((is.na(df_local[[a1]][i]) || df_local[[a1]][i] == "") &&
            !(is.na(df_local[[a2]][i]) || df_local[[a2]][i] == "")) {
            df_local[[a1]][i] <- df_local[[a2]][i]
            df_local[[a2]][i] <- NA
            swapped_total <- swapped_total + 1
            if (!quiet_local) {
              message(
                sprintf(
                  "\t🔄 Reordered %s alleles for sample %s (swapped _1 and _2)",
                  gene,
                  get_sample_id(i)
                )
              )
            }
          }
        }
      }
    }

    if (!quiet_local && swapped_total == 0) {
      message("\tℹ️ No DRB alleles needed reordering.")
    }
    df_local
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Step 2: Fill missing non-DRB alleles
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  fill_other_genes <- function(df_local, quiet_local) {
    allele_cols <- grep("_[12]$", names(df_local), value = TRUE)
    gene_names <- unique(sub("_[12]$", "", allele_cols))
    non_drb_genes <- setdiff(gene_names, c("DRB1", "DRB3", "DRB4", "DRB5"))

    filled_count <- 0
    for (gene in non_drb_genes) {
      a1 <- paste0(gene, "_1")
      a2 <- paste0(gene, "_2")
      if (all(c(a1, a2) %in% names(df_local))) {
        for (i in seq_len(nrow(df_local))) {
          if ((is.na(df_local[[a2]][i]) || df_local[[a2]][i] == "") &&
            !is.na(df_local[[a1]][i]) && df_local[[a1]][i] != "") {
            df_local[[a2]][i] <- df_local[[a1]][i]
            filled_count <- filled_count + 1
            if (!quiet_local) {
              message(
                sprintf(
                  "\t   Filled missing %s_2 allele for sample %s by copying %s_1",
                  gene,
                  get_sample_id(i),
                  gene
                )
              )
            }
          }
        }
      }
    }

    if (!quiet_local && filled_count == 0) {
      message("\tℹ️ No non-DRB alleles needed filling.")
    }
    df_local
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Step 3: DRB1 vs DRB3/4/5 linkage compliance
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  check_and_fill_drb_linkage <- function(df_local, quiet_local) {
    linkage_map <- list(
      "DRB3" = c("03", "11", "12", "13"),
      "DRB4" = c("07", "08", "09"),
      "DRB5" = c("15", "16")
    )

    extract_group <- function(allele) {
      if (is.na(allele) || allele == "") {
        return(NA)
      }
      m <- regmatches(allele, regexpr("\\*(\\d{2})", allele))
      if (length(m) == 0) {
        return(NA)
      }
      sub("\\*", "", m)
    }

    compliance_vec <- logical(nrow(df_local))
    n_noncompliant <- 0

    for (i in seq_len(nrow(df_local))) {
      drb1_1 <- df_local[["DRB1_1"]][i]
      drb1_2 <- df_local[["DRB1_2"]][i]
      drb1_groups <- na.omit(unique(c(extract_group(drb1_1), extract_group(drb1_2))))
      compliant <- TRUE

      for (drb_gene in names(linkage_map)) {
        a1 <- paste0(drb_gene, "_1")
        a2 <- paste0(drb_gene, "_2")
        if (!all(c(a1, a2) %in% names(df_local))) next

        drb1_requires <- any(drb1_groups %in% linkage_map[[drb_gene]])
        alleles_present <- !all(is.na(df_local[i, c(a1, a2)]) | df_local[i, c(a1, a2)] == "")

        if (drb1_requires && !alleles_present) {
          fill_value <- drb1_1
          if (!is.na(fill_value) && fill_value != "") {
            df_local[i, c(a1, a2)] <- fill_value
            if (!quiet_local) {
              message(
                sprintf(
                  "\t   Filled missing %s alleles for sample %s using DRB1_1 allele",
                  drb_gene, get_sample_id(i)
                )
              )
            }
          } else {
            compliant <- FALSE
            if (!quiet_local) {
              message(
                sprintf(
                  "\t⚠️ Sample %s missing required %s alleles and DRB1 allele not available for fill",
                  get_sample_id(i), drb_gene
                )
              )
            }
          }
        } else if (!drb1_requires && alleles_present) {
          compliant <- FALSE
          if (!quiet_local) {
            message(
              sprintf(
                "\t⚠️ Sample %s has %s alleles present but DRB1 group(s) [%s] does not require it",
                get_sample_id(i), drb_gene, paste(drb1_groups, collapse = "/")
              )
            )
          }
        }
      }

      compliance_vec[i] <- compliant
      if (!compliant) n_noncompliant <- n_noncompliant + 1
    }

    if (!quiet_local && n_noncompliant > 0) {
      message(sprintf("\t⚠️ Found %d samples with DRB1-DRBx linkage inconsistencies.", n_noncompliant))
    }

    list(df = df_local, compliance = compliance_vec)
  }

  # Execute the steps
  df <- reorder_drb_alleles(df, quiet)
  df <- fill_other_genes(df, quiet)
  drb_result <- check_and_fill_drb_linkage(df, quiet)

  # Prepare output
  cleaned_df <- drb_result$df
  compliance_flag <- drb_result$compliance

  qc_df <- if (!is.null(id_cols)) {
    cleaned_df[, id_cols, drop = FALSE]
  } else {
    data.frame(RowNumber = seq_len(nrow(cleaned_df)))
  }
  qc_df$DRB145_Compliant <- compliance_flag

  list(cleaned = cleaned_df, qc = qc_df)
}
