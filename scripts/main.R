#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("-m", "--mutations"), type = "character", help = "Path to mutations input file", metavar = "FILE"),
  make_option(c("-s", "--signatures"), type = "character", help = "Path to signatures input file", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory", metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$maf) || is.null(opt$signatures) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("All options must be provided: --maf, --signatures, --output")
}

source("deconstructsigs_assignment.R")
source("sigminer_sigprofiler_assignment.R")

run_deconstructsigs(mutations = opt$mutations, signatures = opt$signatures, output_dir = opt$output)
run_sigminer_sigprofiler(mutations = opt$mutations, signatures = opt$signatures, output_dir = opt$output)

cat("All signature assignments completed successfully.\n")
