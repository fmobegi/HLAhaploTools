#' Reformat Raw HLA Typing Data
#'
#' Cleans and reformats raw HLA typing data for downstream analysis.
#' Includes trimming whitespace, column sorting, and ensuring allele fields
#' are in full gene format (e.g., A*01:01).
#' If it starts with * (e.g., *01:02), prepend the gene name only (e.g., A*01:02).
#' Else if it does not contain * at all (e.g., 01:02), prepend gene + * (e.g., A*01:02).
#' Else (already has gene prefix, e.g., A*01:02), leave as is.
#' Missing second alleles (_2) will be filled with first allele (_1) if present.
#' Rows are reordered by FAMILY_ID and Family_Member with levels F, M, C1, C2 if isFamilyData is TRUE.
#'
#' @param df A data frame or tibble containing raw HLA typing data.
#' @param isFamilyData Logical, default TRUE. If FALSE, Family_Member is ignored and no family-based ordering applied.
#'
#' @return A cleaned and standardized tibble with sorted columns, proper allele format,
#'   filled missing second alleles, and optionally rows ordered by family.
#' @importFrom janitor clean_names
#' @importFrom dplyr across mutate arrange
#' @importFrom forcats fct_relevel
#' @export
reformat_typing_data <- function(df, isFamilyData = TRUE) {
  if (!is.data.frame(df)) {
    stop("Input must be a data frame or tibble.")
  }

  # Preserve original capitalization in column names
  df <- janitor::clean_names(df, case = "none")

  # Trim whitespace in all character columns
  df <- dplyr::mutate(df, dplyr::across(where(is.character), ~ trimws(.)))

  # Define groups
  id_cols  <- c("FAMILY_ID")
  if (isFamilyData) {
    id_cols <- c(id_cols, "Family_Member")
  }

  class1   <- c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2")
  class2   <- c(
    "DRB1_1",
    "DRB1_2",
    "DRB3_1",
    "DRB3_2",
    "DRB4_1",
    "DRB4_2",
    "DRB5_1",
    "DRB5_2",
    "DQA1_1",
    "DQA1_2",
    "DQB1_1",
    "DQB1_2",
    "DPA1_1",
    "DPA1_2",
    "DPB1_1",
    "DPB1_2"
  )
  nonclass <- c("F_1",
                "F_2",
                "G_1",
                "G_2",
                "H_1",
                "H_2",
                "J_1",
                "J_2",
                "E_1",
                "E_2")
  mic      <- c("MICA_1", "MICA_2", "MICB_1", "MICB_2")

  all_allele_cols <- c(class1, class2, nonclass, mic)
  all_order <- c(id_cols, all_allele_cols)
  present   <- intersect(all_order, names(df))
  remaining <- setdiff(names(df), all_order)

  # Ensure allele values include gene name prefix (e.g., A*01:01)
  for (col in intersect(all_allele_cols, names(df))) {
    gene <- sub("_.*", "", col)  # Extract gene from column name

    df[[col]] <- ifelse(
      is.na(df[[col]]) | grepl("^\\s*$", df[[col]]),
      NA_character_,
      ifelse(
        grepl(paste0("^\\*"), df[[col]]),
        # starts with '*'
        paste0(gene, df[[col]]),
        # prepend gene only
        ifelse(
          grepl("\\*", df[[col]]),
          # contains '*' anywhere else
          df[[col]],
          # prepend gene and '*'
          paste0(gene, "*", df[[col]])
        )
      )
    )
  }

  # Fill missing _2 alleles with _1 if _1 is present || Assumes homozygosity
  for (gene in unique(sub("_.*", "", all_allele_cols))) {
    allele1 <- paste0(gene, "_1")
    allele2 <- paste0(gene, "_2")

    if (allele1 %in% names(df) && allele2 %in% names(df)) {
      df[[allele2]] <- ifelse((is.na(df[[allele2]]) |
                                 df[[allele2]] == "") &
                                !is.na(df[[allele1]]) &
                                df[[allele1]] != "", df[[allele1]], df[[allele2]])
    }
  }

  # Reorder columns
  df <- df[, c(present, remaining), drop = FALSE]

  # Reorder rows by FAMILY_ID and Family_Member factor levels if family data
  if (isFamilyData &&
      all(c("Family_Member", "FAMILY_ID") %in% names(df))) {
    desired_order <- c("F", "M", "C1", "C2")
    df <- dplyr::mutate(df,
                        Family_Member = forcats::fct_relevel(Family_Member, desired_order))
    df <- dplyr::arrange(df, FAMILY_ID, Family_Member)
  }

  # Remove completely empty rows and columns
  df <- df[rowSums(is.na(df)) < ncol(df), , drop = FALSE]
  df <- df[, colSums(is.na(df)) < nrow(df), drop = FALSE]

  return(df)
}
