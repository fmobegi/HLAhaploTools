# haplo.em (haplo.stats) on a MAC table, for scripts/bench_haplostats.py.
#
#   Rscript haplostats_em.R input.dat out_dir n_try seed
#
# input.dat: id, then two columns per locus (allele names without the locus
# prefix; an empty cell is an untyped allele). Writes out_dir/hfs.tsv
# (haplotype<TAB>frequency, locus prefixes added) and out_dir/run.tsv
# (lnlike, converge, n.try, subjects used).

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
out <- args[2]
n_try <- as.integer(args[3])
seed <- as.integer(args[4])
lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(lib)) .libPaths(c(path.expand(lib), .libPaths()))
suppressPackageStartupMessages(library(haplo.stats))

header <- strsplit(readLines(input, n = 1), "\t")[[1]]
loci <- header[seq(2, length(header), by = 2)]
tab <- read.delim(input, header = FALSE, skip = 1, colClasses = "character",
                  na.strings = "", quote = "")
geno <- as.matrix(tab[, -1, drop = FALSE])

set.seed(seed)
fit <- haplo.em(geno, locus.label = loci, miss.val = c(NA),
                control = haplo.em.control(n.try = n_try))

haps <- fit$haplotype
names_ <- apply(haps, 1, function(r) paste0(loci, "*", r, collapse = "~"))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
writeLines(paste0(names_, "\t", format(fit$hap.prob, digits = 15, scientific = TRUE)),
           file.path(out, "hfs.tsv"))
writeLines(c(paste0("lnlike\t", format(fit$lnlike, digits = 15)),
             paste0("converge\t", fit$converge),
             paste0("n_try\t", n_try),
             paste0("subjects\t", length(fit$subj.id))),
           file.path(out, "run.tsv"))
