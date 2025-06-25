#' Reformat Raw HLA Typing Data
#'
#' Cleans and standardizes HLA typing data for downstream analysis.
#' This includes trimming whitespace, formatting allele fields (e.g., A*01:01),
#' filling in missing second alleles (_2) using first alleles (_1), and
#' reordering rows by FAMILY_ID and Family_Member (F, M, C1, C2) if desired.
#'
#' @param df A data frame or tibble containing raw HLA typing data.
#' @param isFamilyData Logical. Default TRUE. If FALSE, Family_Member is ignored and no family sorting is applied.
#'
#' @return A tibble with formatted, cleaned, and optionally ordered HLA typing data.
#' @importFrom janitor clean_names
#' @importFrom dplyr across mutate arrange
#' @importFrom forcats fct_relevel
#' @export
reformat_typing_data <- function(df, isFamilyData = TRUE) {
  if (!is.data.frame(df)) {
    stop("Input must be a data frame or tibble.")
  }

  # Clean column names but preserve case
  df <- janitor::clean_names(df, case = "none")
  df <- dplyr::mutate(df, dplyr::across(where(is.character), trimws))

  id_cols <- "FAMILY_ID"
  if (isFamilyData) id_cols <- c(id_cols, "Family_Member")

  class1 <- c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2")
  class2 <- c(
    "DRB1_1", "DRB1_2", "DRB3_1", "DRB3_2", "DRB4_1", "DRB4_2", "DRB5_1", "DRB5_2",
    "DQA1_1", "DQA1_2", "DQB1_1", "DQB1_2", "DPA1_1", "DPA1_2", "DPB1_1", "DPB1_2"
  )
  nonclass <- c("F_1", "F_2", "G_1", "G_2", "H_1", "H_2", "J_1", "J_2", "E_1", "E_2")
  mic <- c("MICA_1", "MICA_2", "MICB_1", "MICB_2")

  all_allele_cols <- c(class1, class2, nonclass, mic)
  all_order <- c(id_cols, all_allele_cols)
  present <- intersect(all_order, names(df))
  remaining <- setdiff(names(df), all_order)

  # Format alleles to include gene name prefix as needed
  for (col in intersect(all_allele_cols, names(df))) {
    gene <- sub("_.*", "", col)

    df[[col]] <- ifelse(
      is.na(df[[col]]) | grepl("^\\s*$", df[[col]]),
      NA_character_,
      ifelse(grepl("^\\*", df[[col]]),
        paste0(gene, df[[col]]), # starts with '*'
        ifelse(grepl("\\*", df[[col]]),
          df[[col]], # already formatted
          paste0(gene, "*", df[[col]]) # bare allele
        )
      )
    )
  }

  # Copy _1 allele to _2 if _2 is blank or missing
  for (gene in unique(sub("_.*", "", all_allele_cols))) {
    a1 <- paste0(gene, "_1")
    a2 <- paste0(gene, "_2")
    if (a1 %in% names(df) && a2 %in% names(df)) {
      df[[a2]] <- ifelse(
        (is.na(df[[a2]]) | df[[a2]] == "") &
          !is.na(df[[a1]]) & df[[a1]] != "",
        df[[a1]],
        df[[a2]]
      )
    }
  }

  # Reorder columns
  df <- df[, c(present, remaining), drop = FALSE]

  # Sort by FAMILY_ID and Family_Member if relevant
  if (isFamilyData && all(c("FAMILY_ID", "Family_Member") %in% names(df))) {
    df <- dplyr::mutate(df,
      Family_Member = forcats::fct_relevel(Family_Member, c("F", "M", "C1", "C2"))
    )
    df <- dplyr::arrange(df, FAMILY_ID, Family_Member)
  }

  # Drop empty rows and columns
  df <- df[rowSums(is.na(df)) < ncol(df), , drop = FALSE]
  df <- df[, colSums(is.na(df)) < nrow(df), drop = FALSE]

  return(df)
}
