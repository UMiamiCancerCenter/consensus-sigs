library(dplyr)
library(maftools)
library(ggplot2)
library(cowplot)
library(reticulate)
library(SigProfilerAssignmentR)
library(reticulate)
library(BSgenome.Hsapiens.UCSC.hg19)
library(quadprog)
library(sigminer)
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-m", "--maf"), type = "character",
              help = "Path to the MAF file"),
  make_option(c("-s", "--signatures"), type = "character",
              help = "Path to the signatures file"),
  make_option(c("-o", "--output"), type = "character",
              help = "Path to the output directory")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if all required options are provided
if (is.null(opt$maf) || is.null(opt$signatures) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("All options must be provided: --maf, --signatures, --output")
}
 
# Get input data
genie_maf<-read_tsv(opt$maf)

results_dir <- file.path(opt$output,"/results")
output_dir <- file.path("./tmp/output")

dir.create(file.path(results_dir, "SigProfiler_Cosmic_V3.4"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "SigMiner_Cosmic_V3.4"), recursive = TRUE, showWarnings = FALSE)


# Sigminer has a great tool for loading and extracting the Tally from maf data based on maftool package
MAF <- sigminer::read_maf(maf = opt$maf)

#Tally the SNV in 96 possible alterations   
mt_tally_SNV <- sigminer::sig_tally(MAF,
                                    ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                                    genome_build = "hg19",
                                    useSyn = TRUE)

saveRDS(MAF, file.path(output_dir, "RObjects", "MAF_all.rds"))

saveRDS(mt_tally_SNV, file.path(output_dir, "RObjects",  "mt_tally_SNV.rds"))

# Save a transposed version for SigProfiler
transposed <- data.frame(MutationType = colnames(mt_tally_SNV$nmf_matrix), t(mt_tally_SNV$nmf_matrix))

data.table::fwrite(x = transposed, sep = "\t", quote = F, row.names = F, file = file.path(output_dir, "ALL_nmf_matrix.txt"))

data <- file.path(output_dir, "ALL_nmf_matrix.txt")

cosmic_fit(data, 
           output = file.path(results_dir, "SigProfiler_Cosmic_V3.4"), 
           input_type='matrix', 
           context_type="96",
           collapse_to_SBS96=TRUE, 
           cosmic_version=3.4, 
           exome=TRUE,
           genome_build="GRCh37", 
           signature_database=NULL,
           exclude_signature_subgroups=NULL, 
           export_probabilities=TRUE,
           export_probabilities_per_mutation=FALSE, 
           make_plots=TRUE,
           sample_reconstruction_plots=FALSE, 
           verbose=T, 
           nnls_add_penalty   = 0.05,
           nnls_remove_penalty= 0.0,    
           initial_remove_penalty=0.0)


COSMIC_V3 <- sigminer::sig_fit(catalogue_matrix = t(mt_tally_SNV$nmf_matrix), 
                               sig_db = "latest_SBS_GRCh37", sig_index = "ALL", 
                               return_class = "data.table", exome=TRUE)

data.table::fwrite(sep = ",", x = COSMIC_V3, file = file.path(results_dir, "SigMiner_Cosmic_V3.4", "SigMiner.csv"))
