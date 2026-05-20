#' HLAhaploTools: Tools for Inference and Analysis of Extended HLA Haplotypes
#'
#' Provides tools for constructing, analysing, and visualising extended
#' human leukocyte antigen (HLA) haplotypes, including family‑based
#' segregation analysis, EM‑based inference, and population‑level summaries.
#'
#' @name HLAhaploTools
#' @aliases run_haplotype_analysis
#' @keywords package
#'
#' @importFrom scales percent label_percent alpha
#' @importFrom stats reorder frequency na.omit
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate distinct slice_max relocate select filter
#' @importFrom tibble tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_y_continuous
#' @importFrom haplo.stats haplo.em
#' @import dplyr
#' @import forcats
#' @import furrr
#' @import future
#' @import future.apply
#' @import ggplot2
#' @import HLAtools
#' @import immunotation
#' @import janitor
#' @import purrr
#' @import RColorBrewer
#' @import readr
#' @import readxl
#' @import tibble
#' @import tidyr
#' @import patchwork
"_PACKAGE"
