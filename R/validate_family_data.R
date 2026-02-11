# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# validate_family_data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Validate Family Data Structure
#'
#' Checks whether the input dataset reflects a valid family study format and
#' provides structural insights. This includes detecting father, mother, and child roles
#' based on the `Family_Member` column and reporting completeness of family trios.
#'
#' @param df A data frame or tibble containing HLA typing data.
#' @param stop_if_invalid Logical. If TRUE, stops execution if data is not valid family format.
#' @param verbose Logical. If TRUE, prints detailed validation output.
#'
#' @return Logical `TRUE` if valid family structure is detected, `FALSE` otherwise.
#' @keywords internal
validate_family_data <- function(df,
                                 stop_if_invalid = TRUE,
                                 verbose = TRUE) {
   # Detect data type using quiet mode
   detect <- detect_data_type(df, quiet = !verbose)
   is_family <- detect$is_family

   if (!is_family) {
      if (stop_if_invalid) {
         stop(
            "\n\t❌ Error: Data does not appear to be in family study format.\n",
            "\t   Use other functions for registry or population data."
         )
      }
      return(FALSE)
   }

   if (verbose) {
      message("\n\t🔍 Validating family structure...")

      family_col <- "FAMILY_ID"
      member_col <- "Family_Member"

      # Summary counts
      n_families <- length(unique(df[[family_col]]))
      n_fathers <- sum(df[[member_col]] %in% c("F", "Father"), na.rm = TRUE)
      n_mothers <- sum(df[[member_col]] %in% c("M", "Mother"), na.rm = TRUE)
      n_children <- sum(grepl("^C\\d+$|Child", df[[member_col]]), na.rm = TRUE)

      message(sprintf(
         "\t   Found %d families with %d fathers, %d mothers, %d children.",
         n_families, n_fathers, n_mothers, n_children
      ))

      # Trio completeness check
      complete_families <- sum(sapply(unique(df[[family_col]]), function(fam) {
         fam_data <- df[df[[family_col]] == fam, ]
         has_father <- any(fam_data[[member_col]] %in% c("F", "Father"))
         has_mother <- any(fam_data[[member_col]] %in% c("M", "Mother"))
         has_child <- any(grepl("^C\\d+$|Child", fam_data[[member_col]]))
         has_father && has_mother && has_child
      }))

      message(sprintf(
         "\t   %d/%d families have complete trios (father, mother, child).",
         complete_families, n_families
      ))
   }

   TRUE
}
