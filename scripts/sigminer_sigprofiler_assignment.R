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

python_path <- Sys.getenv("RETICULATE_PYTHON", unset = NA)
use_python(python_path, required = TRUE)

run_sigminer_sigprofiler <- function(maf, signatures, results_dir) {
  tmp_dir <- file.path("./tmp/output")

  dir.create(file.path(tmp_dir, "RObjects"), recursive = TRUE, showWarnings = FALSE)

  maf_table <- sigminer::read_maf(maf = maf)

  # Tally the SNV in 96 possible alterations
  mt_tally_snv <- sigminer::sig_tally(maf_table,
    ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
    genome_build = "hg19",
    useSyn = TRUE
  )

  saveRDS(maf_table, file.path(tmp_dir, "RObjects", "MAF_all.rds"))
  saveRDS(mt_tally_snv, file.path(tmp_dir, "RObjects", "mt_tally_SNV.rds"))

  # Save a transposed version for SigProfiler
  transposed <- data.frame(MutationType = colnames(mt_tally_snv$nmf_matrix), t(mt_tally_snv$nmf_matrix))

  data.table::fwrite(x = transposed, sep = "\t", quote = FALSE, row.names = FALSE, file = file.path(tmp_dir, "ALL_nmf_matrix.txt"))

  data <- file.path(tmp_dir, "ALL_nmf_matrix.txt")

  # Run SigProfiler assignment
  cosmic_fit(data,
    output = file.path(results_dir, "SigProfiler_Cosmic_V3.4"),
    input_type = "matrix",
    context_type = "96",
    collapse_to_SBS96 = TRUE,
    cosmic_version = 3.4,
    exome = TRUE,
    genome_build = "GRCh37",
    signature_database = NULL,
    exclude_signature_subgroups = NULL,
    export_probabilities = TRUE,
    export_probabilities_per_mutation = FALSE,
    make_plots = TRUE,
    sample_reconstruction_plots = FALSE,
    verbose = TRUE,
    nnls_add_penalty = 0.05,
    nnls_remove_penalty = 0.0,
    initial_remove_penalty = 0.0
  )

  # Run SigMiner assignment
  cosmic_v3 <- sigminer::sig_fit(
    catalogue_matrix = t(mt_tally_snv$nmf_matrix),
    sig_db = "latest_SBS_GRCh37", sig_index = "ALL",
    return_class = "data.table", exome = TRUE
  )

  data.table::fwrite(sep = ",", x = cosmic_v3, file = file.path(results_dir, "SigMiner_Cosmic_V3.4", "SigMiner.csv"))
}
