library(readr)
library(deconstructSigs)
library(dplyr)
library(tibble)
library(ggplot2)
library(forcats)
library(rlang)
 
run_deconstructsigs <- function(maf, signatures, result_dir) {
  # Get input data
  maf_table<-read_tsv(maf)
  cosmic_v3_4_signatures<-readRDS(signatures)

  # Assign signatures
  sample_ids<-unique(maf_table$Tumor_Sample_Barcode)
  sample_sigs<-list()
  failed<-c()
  n<- 0
  sample_signatures<-lapply(sample_ids, function(sample_id){
    tryCatch({
      one_sample_mutations<-maf_table%>%
        filter(
          .data$Tumor_Sample_Barcode == sample_id
        )%>%
        select(
          Sample = "Tumor_Sample_Barcode", chr = "Chromosome", pos = "Start_Position", ref = "Reference_Allele", alt = "Tumor_Seq_Allele2"
        )
      
      sigs_input<-mut.to.sigs.input(mut.ref = as.data.frame(one_sample_mutations),
                                    sample.id = "Sample",
                                    chr = "chr", 
                                    pos = "pos", 
                                    ref = "ref", 
                                    alt = "alt")
      
      result <- whichSignatures(tumor.ref = sigs_input, signatures.ref = cosmic_v3_4_signatures, tri.counts.method = 'exome', contexts.needed = TRUE)
      return(result)
      }, error = function(e) {
        failed<-c(failed, sample_id)
        return(NULL)
      })
    })
  names(sample_signatures)<-sample_ids

  # Get weights
  signature_weights_list<-lapply(names(sample_signatures), function(sample_id){
    weights<-sample_signatures[[sample_id]]$weights
    return(weights)
  })
  names(signature_weights_list)<-names(sample_signatures)


  # All signatures matrix
  all_weights_matrix<-as.matrix(do.call(rbind, signature_weights_list))
  write_csv(rownames_to_column(as.data.frame(all_weights_matrix), var = "sample_name"), file.path(result_dir, "deconstructSigs_V3.4", "deconstructSigs_cosmic_v3_4_all_signatures.csv"))
}
