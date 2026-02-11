# ~~~~~~~~~~~~~~~~~~~~~
# check_deleted_alleles
# ~~~~~~~~~~~~~~~~~~~~~
#' Check HLA dataframe for deleted alleles
#'
#' @description Analyzes a dataframe containing HLA allele designations, checking for
#' deleted alleles according to the IMGT/HLA database. Downloads and parses the official
#' deleted allele list, then matches any occurrences in the user's data.
#'
#' @param df_hla A dataframe containing HLA allele typing data.
#' @param deleted_file_path Path or URL to the Deleted_alleles.txt file
#'        (default: IMGT/HLA database URL).
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A dataframe containing information about deleted alleles found in the input.
#'
#' @examples
#' \dontrun{
#' results <- check_deleted_alleles(df_hla)
#' }
#' @keywords internal
#' @importFrom curl curl
check_deleted_alleles <- function(
  df_hla,
  deleted_file_path = "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Deleted_alleles.txt",
  quiet = FALSE
) {
   if (!quiet) {
      message("\t🔍 Checking for deleted HLA alleles...")
      message("\t📥 Downloading current list of deleted alleles from:\n\t ", deleted_file_path)
   }

   deleted_alleles <- list()

   tryCatch(
      {
         con <- if (grepl("^https?://", deleted_file_path)) {
            curl::curl(deleted_file_path)
         } else {
            file(deleted_file_path, open = "r")
         }

         lines <- readLines(con, warn = FALSE)
         close(con)

         header_found <- FALSE

         for (line in lines) {
            if (grepl("^#", line) || nchar(trimws(line)) == 0) next
            if (grepl("^AlleleID,Allele,Description", line)) {
               header_found <- TRUE
               next
            }
            if (header_found) {
               split_line <- strsplit(line, ",")[[1]]
               if (length(split_line) >= 2) {
                  allele_name <- trimws(split_line[2])
                  description <- if (length(split_line) >= 3) {
                     paste(split_line[3:length(split_line)], collapse = ",")
                  } else {
                     "No description"
                  }
                  deleted_alleles[[allele_name]] <- description
               }
            }
         }

         if (!quiet) {
            message(paste0("\t📊 Loaded ", length(deleted_alleles), " deleted alleles."))
         }
      },
      error = function(e) {
         message("\t❌ Error reading deleted alleles file: ", e$message)
         data.frame()
      }
   )

   if (length(deleted_alleles) == 0) {
      if (!quiet) {
         message("\t⚠️ No deleted alleles found in the reference file.")
      }
      data.frame()
   }

   results <- data.frame(
      FAMILY_ID = character(),
      Family_Member = character(),
      Column = character(),
      Deleted_Allele = character(),
      Description = character(),
      stringsAsFactors = FALSE
   )

   hla_columns <- setdiff(names(df_hla), c("FAMILY_ID", "Family_Member"))

   for (i in seq_len(nrow(df_hla))) {
      for (col in hla_columns) {
         allele <- df_hla[i, col]
         if (is.na(allele) || allele == "") next
         if (allele %in% names(deleted_alleles)) {
            results <- rbind(results, data.frame(
               FAMILY_ID = df_hla$FAMILY_ID[i],
               Family_Member = df_hla$Family_Member[i],
               Column = col,
               Deleted_Allele = allele,
               Description = deleted_alleles[[allele]],
               stringsAsFactors = FALSE
            ))

            if (!quiet) {
               message(sprintf(
                  "\t🔴 Found deleted allele: %s in FAMILY_ID=%s, Family_Member=%s, Column=%s",
                  allele, df_hla$FAMILY_ID[i], df_hla$Family_Member[i], col
               ))
            }
         }
      }
   }

   if (!quiet) {
      if (nrow(results) == 0) {
         message("\t✅ No deleted alleles found in the data.")
      } else {
         message(sprintf("\t⚠️ Found %d instances of deleted alleles.", nrow(results)))
      }
   }

   results
}
