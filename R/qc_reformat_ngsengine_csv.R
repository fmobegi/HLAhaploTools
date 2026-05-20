# path <- "W:/Pathology/FSH/Immunology/Validation Data and Research Projects (AD13)/19 IHWS projects/Full gene HLA haplotypes/Data files"
#
# all_files <- list.files(
#    path,
#    pattern = "\\.csv$",
#    full.names = TRUE
# )
#
# process_one <- function(f) {
#    readr::read_csv(f, show_col_types = FALSE) %>%
#       dplyr::select(SampleID, Locus, NMDP1, NMDP2) %>%
#       tidyr::pivot_wider(
#          id_cols = SampleID,
#          names_from = Locus,
#          values_from = c(NMDP1, NMDP2),
#          names_glue = "{Locus}_{.value}"
#       ) %>%
#       dplyr::rename_with(~ sub("NMDP1", "1", .x), contains("NMDP1")) %>%
#       dplyr::rename_with(~ sub("NMDP2", "2", .x), contains("NMDP2"))
# }
#
# gene_order <- c(
#    "A_1", "A_2", "B_1", "B_2", "C_1", "C_2",
#    "E_1", "E_2", "F_1", "F_2", "G_1", "G_2", "H_1", "H_2", "K_1", "K_2",
#    "MICA_1", "MICA_2", "MICB_1", "MICB_2",
#    "DRB1_1", "DRB1_2", "DRB3_1", "DRB3_2", "DRB4_1", "DRB4_2", "DRB5_1", "DRB5_2",
#    "DQA1_1", "DQA1_2", "DQB1_1", "DQB1_2", "DPA1_1", "DPA1_2", "DPB1_1", "DPB1_2"
# )
#
# final_df <- purrr::map_dfr(all_files, process_one) %>%
#    dplyr::select(SampleID, any_of(gene_order), everything()) %>%
#    fill_missing_alleles()
#
# outfile <- paste0(path, "/Combined_results_tabular.csv.txt")
# write_tsv(final_df$cleaned, outfile)
