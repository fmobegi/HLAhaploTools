# ~~~~~~~~~~~~~~~~
# trim_hla_results
# ~~~~~~~~~~~~~~~~
#' Trim Resolution of HLA Alleles Across a Data Frame
#'
#' Applies `HLAtools::multiAlleleTrim()` to all columns in a tibble or data frame,
#' reducing allele resolution (e.g., from 4-field to 2-field notation).
#'
#' Handles MAC-decoded alleles containing "/" separators by trimming each component.
#' NA values and empty strings are preserved, and trimming failures are gracefully handled with fallback warnings.
#'
#' @param df A tibble or data frame with HLA alleles as character columns.
#' @param resolution Integer. Desired resolution level (e.g., 2 or 4). Default is 3.
#' @param append Logical. If TRUE (default), ambiguity suffixes (e.g., "G", "N") are retained after trimming.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A tibble with trimmed allele values.
#'
#' @examples
#' \dontrun{
#' # Load and trim alleles to 2-field resolution
#' typed_data <- load_typing_data("hla_typing.csv") %>%
#'    reformat_typing_data() %>%
#'    trim_hla_results(resolution = 2)
#' }
#' @importFrom HLAtools multiAlleleTrim
#' @importFrom tibble as_tibble
#'
#' @export
trim_hla_results <- function(df,
                             resolution = 3,
                             append = TRUE,
                             quiet = FALSE) {
   # Check for HLAtools package
   if (!requireNamespace("HLAtools", quietly = TRUE)) {
      stop(
         "\t❌ The 'HLAtools' package is required. Please install it from GitHub: remotes::install_github('DKMS/HLAtools')"
      )
   }

   if (!quiet) {
      message(sprintf("\t✂️ Trimming alleles to %d-field resolution...", resolution))
   }

   # Count of trimmed alleles for reporting
   trim_count <- 0
   fail_count <- 0
   mac_count <- 0

   # Apply trimming to each column
   df[] <- lapply(df, function(col) {
      if (!is.character(col)) {
         return(col)
      }

      sapply(col, function(allele) {
         if (is.na(allele) || allele == "") {
            return(allele)
         }

         # Check if this is a MAC-decoded allele containing "/"
         if (grepl("/", allele)) {
            mac_count <<- mac_count + 1

            # Split by "/" and trim each component
            components <- unlist(strsplit(allele, "/"))
            trimmed_components <- sapply(components, function(component) {
               tryCatch(
                  {
                     trimmed <- HLAtools::multiAlleleTrim(component,
                        resolution = resolution,
                        append = append
                     )
                     trim_count <<- trim_count + 1
                     return(trimmed)
                  },
                  error = function(e) {
                     fail_count <<- fail_count + 1
                     if (!quiet) {
                        message(
                           sprintf(
                              "\t⚠️ Trimming failed for '%s' in MAC-decoded string — keeping original.",
                              component
                           )
                        )
                     }
                     return(component)
                  }
               )
            })

            # Recombine with "/"
            return(paste(trimmed_components, collapse = "/"))
         }

         # Regular allele (no MAC decoding)
         tryCatch(
            {
               trimmed <- HLAtools::multiAlleleTrim(allele, resolution = resolution, append = append)
               trim_count <<- trim_count + 1
               return(trimmed)
            },
            error = function(e) {
               fail_count <<- fail_count + 1
               if (!quiet) {
                  message(sprintf(
                     "\t⚠️ Trimming failed for '%s' — keeping original.",
                     allele
                  ))
               }
               return(allele)
            }
         )
      })
   })

   if (!quiet) {
      message(
         sprintf(
            "\t✅ Trimmed %d alleles to %d-field resolution.",
            trim_count - fail_count,
            resolution
         )
      )
      if (mac_count > 0) {
         message(
            sprintf(
               "\t🧩 Processed %d MAC-decoded alleles with multiple options.",
               mac_count
            )
         )
      }
      if (fail_count > 0) {
         message(sprintf("\t⚠️ %d alleles could not be trimmed properly.", fail_count))
      }
   }
   df <- tibble::as_tibble(df)
   df
}
