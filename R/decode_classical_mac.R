# ~~~~~~~~~~~~~~~~~~~~
# decode_classical_mac
# ~~~~~~~~~~~~~~~~~~~~
#' Decode MAC-Encoded Alleles in Classical HLA Columns
#'
#' Uses immunotation to decode multi-allele codes (MAC) in classical HLA loci.
#'
#' @param df A data frame with HLA allele columns (e.g. A_1, B_2, etc).
#' @param loci Optional set of classical genes to decode.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return Data frame with decoded MAC values.
#' @examples
#' \dontrun{
#' # Decode all classical loci
#' typed_data <- load_typing_data("hla_typing.csv") %>%
#'   reformat_typing_data() %>%
#'   decode_classical_mac()
#'
#' # Only decode Class I loci
#' class1_decoded <- decode_classical_mac(typed_data, loci = c("A", "B", "C"))
#' }
#' @export
#' @importFrom immunotation decode_MAC
decode_classical_mac <- function(df,
                                 loci = c(
                                   "A",
                                   "B",
                                   "C",
                                   "DRB1",
                                   "DRB3",
                                   "DRB4",
                                   "DRB5",
                                   "DQA1",
                                   "DQB1",
                                   "DPA1",
                                   "DPB1"
                                 ),
                                 quiet = FALSE) {
  # Check for immunotation package
  if (!requireNamespace("immunotation", quietly = TRUE)) {
    stop(
      "\t❌ The 'immunotation' package is required. Please install it with BiocManager::install('immunotation')."
    )
  }

  if (!quiet) {
    message("\t🔍 Decoding MAC strings in HLA alleles...")
  }

  # Find columns to process
  allele_cols <- intersect(unlist(lapply(loci, function(l) {
    paste0(l, c(
      "_1", "_2"
    ))
  })), names(df))

  if (length(allele_cols) == 0) {
    warning("\t⚠️ No matching HLA columns found for MAC decoding.")
    return(df)
  }

  # MAC pattern
  mac_pattern <- "^[A-Z0-9]+\\*[0-9]{2}:[A-Z]+$"
  mac_count <- 0

  # Process each column
  for (col in allele_cols) {
    df[[col]] <- sapply(df[[col]], function(allele) {
      if (is.na(allele) || allele == "") {
        return(allele)
      }
      if (grepl(mac_pattern, allele)) {
        mac_count <<- mac_count + 1
        tryCatch(
          immunotation::decode_MAC(allele),
          error = function(e) {
            warning(sprintf("\t⚠️ Failed to decode MAC: %s", allele))
            allele
          }
        )
      } else {
        allele
      }
    })
  }

  if (!quiet) {
    if (mac_count > 0) {
      message(sprintf("\t✅ Decoded %d MAC strings.", mac_count))
    } else {
      message("\tℹ️ No MAC strings found to decode.")
    }
  }

  df
}
