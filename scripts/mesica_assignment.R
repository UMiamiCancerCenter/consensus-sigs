library(dplyr)
library(tibble)
library(ggplot2)
library(reshape2)
library(rlang)
source("./utils/mesica_utils.R")


run_mesica <- function(maf, results_dir) {
  # Get TCGA embeddings
  fd <- read.csv("./resources/mesica/sig_label_feature_dict.csv") %>% dplyr::arrange(X0)
  pm <- read.csv("./resources/mesica/sig_label_patient_mapping.csv")
  p2c <- read.csv("./resources/mesica/sig_label_pattern2col.csv") %>% dplyr::arrange(X0)
  embeddings <- read.csv("./resources/mesica/sig_label_embeddings.csv")[, -1]
  rownames(embeddings) <- fd$X
  rownames(embeddings)[1:1536] <- p2c$X

  maf_table <- read.table(maf, sep = "\t", header = TRUE, quote = "", check.names = FALSE)

  # Filter redundant samples
  maf_mut <- create_muts(maf_table %>% filter(.data$Reference_Allele != .data$Tumor_Seq_Allele2))

  custom_sig_names <- c("APOBEC", "Clock_SBS1", "Clock_SBS5", "HRD", "Tobacco", "MMR", "UV", "POLE", "other")

  preds_agg_genie_all_list <- cos_pred_fun_efficient(
    sig_names = custom_sig_names,
    emb = embeddings,
    unique(maf_mut$Tumor_Sample_Barcode),
    maf_mut,
    creating_muts = FALSE,
    vep = TRUE,
    emb_len = ncol(embeddings)
  )

  preds_agg_genie_all <- preds_agg_genie_all_list[[3]]
  data.table::fwrite(sep = ",", x = preds_agg_genie_all, file = file.path(results_dir, "MESiCA_Cosmic_V3.4", "MESiCA_assignment.csv"))
}
