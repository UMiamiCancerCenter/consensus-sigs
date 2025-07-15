#!/usr/bin/env Rscript

library(optparse)
source("deconstructsigs_assignment.R")
source("sigminer_sigprofiler_assignment.R")
source("mesica_assignment.R")
source("consensus_pipeline.R")

option_list <- list(
  make_option(c("-m", "--maf"), type = "character", help = "Path to MAF input file", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory", metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$maf) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("All options must be provided: --maf, --output")
}

results_dir <- file.path(opt$output, "consensus_sigs_results")
for (subdir in c("SigProfiler_Cosmic_V3.4", "SigMiner_Cosmic_V3.4", "deconstructSigs_Cosmic_V3.4", "MESiCA_Cosmic_V3.4", "processed", "consensus")) {
  dir.create(file.path(results_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}

run_deconstructsigs(maf = opt$maf, results_dir = results_dir)
run_sigminer_sigprofiler(maf = opt$maf, results_dir = results_dir)
run_mesica(maf = opt$maf, results_dir = results_dir)
run_consensus_pipeline(results_dir = results_dir)

cat("All signature assignments completed successfully.\n")
