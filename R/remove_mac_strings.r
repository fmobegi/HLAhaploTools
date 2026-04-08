# ~~~~~~~~~~~~~~~~
# remove mac strings
# ~~~~~~~~~~~~~~~~
#' Trim Resolution of HLA Alleles Across a Data Frame
#'
#' Applies `HLAtools::multiAlleleTrim()` to all columns in a tibble or data frame,
#' reducing allele resolution (e.g., from 4-field to 2-field notation).
#'
#' Handles MAC-decoded alleles containing "/" separators by trimming each component.
#' Return only first allele.
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
#'   reformat_typing_data() %>%
#'   remove_mac_strings(resolution = 2)
#' }
#' @importFrom HLAtools multiAlleleTrim
#' @importFrom tibble as_tibble
#'
#' @export
remove_mac_strings <- function(df,
                               resolution = 3,
                               append = TRUE,
                               quiet = FALSE) {
  if (!requireNamespace("HLAtools", quietly = TRUE)) {
    stop(
      "\t❌ The 'HLAtools' package is required. Please install it from GitHub: remotes::install_github('DKMS/HLAtools')"
    )
  }

  if (!quiet) {
    message(sprintf("\t✂️ Trimming alleles to %d-field resolution...", resolution))
  }

  trim_count <- 0
  fail_count <- 0
  mac_count <- 0

  df[] <- lapply(df, function(col) {
    if (!is.character(col)) {
      return(col)
    }

    sapply(col, function(allele) {
      if (is.na(allele) || allele == "") {
        return(allele)
      }

      if (grepl("/", allele)) {
        mac_count <<- mac_count + 1

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

        return(trimmed_components[1])
      }

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
