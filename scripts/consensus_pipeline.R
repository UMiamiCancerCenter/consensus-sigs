library(dplyr)
library(rlang)
library(purrr)
source("./utils/consensus_pipeline_utils.R")

# Load algorithm assignments
run_consensus_pipeline <- function(results_dir) {
  decontructsigs_sigs <- data.table::fread(
    sep = ",", header = TRUE,
    file = file.path(results_dir, "deconstructSigs_Cosmic_V3.4", "deconstructSigs_assignment.csv")
  ) %>%
    data.table::setnames("sample_name", "SAMPLE_ID") %>%
    dplyr::arrange(SAMPLE_ID)

  sigprofiler_sigs <- data.table::fread(
    sep = "\t", header = TRUE,
    file = file.path(results_dir, "SigProfiler_Cosmic_V3.4/Assignment_Solution/Activities/Assignment_Solution_Activities.txt")
  ) %>%
    data.table::setnames("Samples", "SAMPLE_ID") %>%
    dplyr::mutate(SAMPLE_ID = gsub("\\.", "-", SAMPLE_ID)) %>%
    dplyr::arrange(SAMPLE_ID)

  sigminer_sigs <- data.table::fread(
    sep = ",", header = TRUE,
    file = file.path(results_dir, "SigMiner_Cosmic_V3.4", "SigMiner_assignment.csv")
  ) %>%
    data.table::setnames("sample", "SAMPLE_ID") %>%
    dplyr::arrange(SAMPLE_ID)
  sigminer_sigs <- sigminer_sigs[rowSums(sigminer_sigs[, -1]) > 0, ]

  mesica_sigs <- data.table::fread(
    sep = ",", header = TRUE,
    file = file.path(results_dir, "MESiCA_Cosmic_V3.4", "MESiCA_assignment.csv")
  ) %>%
    data.table::setnames("sample", "SAMPLE_ID") %>%
    dplyr::arrange(SAMPLE_ID) %>%
    dplyr::select(-other)

  # Binarize the signature assignments
  bin_decontructsigs_sigs <- decontructsigs_sigs %>%
    mutate(across(-any_of("SAMPLE_ID"), ~ as.integer(. > 0))) %>%
    as.data.frame()

  bin_sigminer_sigs <- sigminer_sigs %>%
    mutate(across(-any_of("SAMPLE_ID"), ~ as.integer(. > 0))) %>%
    as.data.frame()

  bin_sigprofiler_sigs <- sigprofiler_sigs %>%
    mutate(across(-any_of("SAMPLE_ID"), ~ as.integer(. > 0))) %>%
    as.data.frame()

  bin_mesica_sigs <- mesica_sigs %>%
    mutate(across(-any_of("SAMPLE_ID"), ~ ifelse(. == "Yes", 1L, ifelse(. == "No", 0L, NA_integer_)))) %>%
    as.data.frame()

  # Select only the top signature per sample and binarize.

  top_decontructsigs_sigs <- data.frame(SAMPLE_ID = decontructsigs_sigs$SAMPLE_ID, ifelse(decontructsigs_sigs[, -1] == apply(decontructsigs_sigs[, -1], 1, max), 1, 0))

  top_sigminer_sigs <- data.frame(SAMPLE_ID = sigminer_sigs$SAMPLE_ID, ifelse(sigminer_sigs[, -1] == apply(sigminer_sigs[, -1], 1, max), 1, 0))

  top_sigprofiler_sigs <- data.frame(SAMPLE_ID = sigprofiler_sigs$SAMPLE_ID, ifelse(sigprofiler_sigs[, -1] == apply(sigprofiler_sigs[, -1], 1, max), 1, 0))

  data.table::fwrite(x = bin_decontructsigs_sigs, sep = ",", row.names = FALSE, file = file.path(results_dir, "processed", "binarized_deconstructSigs.csv"))
  data.table::fwrite(x = bin_sigminer_sigs, sep = ",", row.names = FALSE, file = file.path(results_dir, "processed", "binarized_SigMiner.csv"))
  data.table::fwrite(x = bin_sigprofiler_sigs, sep = ",", row.names = FALSE, file = file.path(results_dir, "processed", "binarized_SigProfiler.csv"))

  data.table::fwrite(x = top_decontructsigs_sigs, sep = ",", row.names = FALSE, file = file.path(results_dir, "processed", "top_Binarized_deconstructSigs.csv"))
  data.table::fwrite(x = top_sigminer_sigs, sep = ",", row.names = FALSE, file = file.path(results_dir, "processed", "top_Binarized_SigMiner.csv"))
  data.table::fwrite(x = top_sigprofiler_sigs, sep = ",", row.names = FALSE, file = file.path(results_dir, "processed", "top_Binarized_SigProfiler.csv"))

  decontructsigs_labeled <- create_mesica_table(bin_decontructsigs_sigs, subgroup_to_sigs)
  decontructsigs_labeled_top <- create_mesica_table(top_decontructsigs_sigs, subgroup_to_sigs)

  sigminer_labeled <- create_mesica_table(bin_sigminer_sigs, subgroup_to_sigs)
  sigminer_labeled_top <- create_mesica_table(top_sigminer_sigs, subgroup_to_sigs)

  sigprofiler_labeled <- create_mesica_table(bin_sigprofiler_sigs, subgroup_to_sigs)
  sigprofiler_labeled_top <- create_mesica_table(top_sigprofiler_sigs, subgroup_to_sigs)

  data.table::fwrite(x = decontructsigs_labeled, sep = ",", row.names = FALSE, file = file.path(results_dir, "processed", "deconstructSigs_labeled.csv"))
  data.table::fwrite(x = decontructsigs_labeled_top, sep = ",", row.names = FALSE, file = file.path(results_dir, "processed", "deconstructSigs_labeled_top.csv"))

  data.table::fwrite(x = sigminer_labeled, sep = ",", row.names = FALSE, file = file.path(results_dir, "processed", "SigMiner_labeled.csv"))
  data.table::fwrite(x = sigminer_labeled_top, sep = ",", row.names = FALSE, file = file.path(results_dir, "processed", "SigMiner_labeled_top.csv"))

  data.table::fwrite(x = sigprofiler_labeled, sep = ",", row.names = FALSE, file = file.path(results_dir, "processed", "SigProfiler_labeled.csv"))
  data.table::fwrite(x = sigprofiler_labeled_top, sep = ",", row.names = FALSE, file = file.path(results_dir, "processed", "SigProfiler_labeled_top.csv"))

  # Consensus Table
  assignment_tables <- list(decontructsigs_labeled_top, sigminer_labeled_top, sigprofiler_labeled_top, bin_mesica_sigs)

  all_samples <- assignment_tables %>%
    map(~ pull(.x, SAMPLE_ID)) %>%
    reduce(union) %>%
    sort()

  assignment_tables_padded <- assignment_tables %>%
    map(
      ~ data.frame(SAMPLE_ID = all_samples) %>%
        left_join(.x, by = "SAMPLE_ID") %>%
        mutate(across(-SAMPLE_ID, ~ ifelse(is.na(.x), 0, .x)))
    )

  assignment_tables_padded <- assignment_tables_padded %>%
    map(
      ~ {
        rownames(.x) <- .x$SAMPLE_ID
        .x <- .x[, -1]
        as.matrix(.x)
      }
    )
  consensus_table <- Reduce(`+`, assignment_tables_padded) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("SAMPLE_ID")

  bin_consensus_table <- consensus_table %>%
    mutate(across(-any_of("SAMPLE_ID"), ~ if_else(.x >= 2, 1, 0)))

  data.table::fwrite(x = consensus_table, sep = ",", row.names = FALSE, file = file.path(results_dir, "consensus", "consensus_table.csv"))
  data.table::fwrite(x = bin_consensus_table, sep = ",", row.names = FALSE, file = file.path(results_dir, "consensus", "bin_consensus_table.csv"))
}
