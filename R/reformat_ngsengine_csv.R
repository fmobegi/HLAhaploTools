library(tidyverse)

path <- "C:/Users/he144508/Downloads/"

all_files <- list.files(
   path,
   pattern = "\\.csv$",
   full.names = TRUE
)

process_one <- function(f) {
   readr::read_csv(f, show_col_types = FALSE) %>%
      dplyr::select(SampleID, Locus, NMDP1, NMDP2) %>%
      tidyr::pivot_wider(
         id_cols = SampleID,
         names_from = Locus,
         values_from = c(NMDP1, NMDP2),
         names_glue = "{Locus}_{.value}"
      ) %>%
      dplyr::rename_with(~ sub("NMDP1", "1", .x), contains("NMDP1")) %>%
      dplyr::rename_with(~ sub("NMDP2", "2", .x), contains("NMDP2"))
}

final_df <- purrr::map_dfr(all_files, process_one)
