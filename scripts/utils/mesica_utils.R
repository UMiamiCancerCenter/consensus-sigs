library(rtracklayer)
library(parallel)
library(reshape2)
library(ggplot2)
library(Biostrings)
library(plyr)
library(Rsamtools)
library(VariantAnnotation)
library(dplyr)
library(survminer)
library(survival)
library(foreach)
library(doParallel)

## Nucleotide Context extraction
create_muts <- function(mut){
  mut_snv <- mut %>% filter(Variant_Type == "SNP")
  chromosomes_snv <-  paste0("chr",mut_snv$Chromosome)
  vr <- VRanges(seqnames = chromosomes_snv,
                ranges = IRanges(start = mut_snv$Start_Position-3, end = mut_snv$End_Position+3),
                ref = mut_snv$Reference_Allele, alt = mut_snv$Tumor_Seq_Allele2)
  t<-SomaticSignatures::mutationContext(vr, BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, k = 7)
  
  mut_snv$alteration <- as.character(t$alteration)
  mut_snv$context <- as.character(t$context)
  mut_snv$upstream <- substr(mut_snv$context,1,3)
  mut_snv$downstream <- substr(mut_snv$context,5,7)
  mut_snv$var_type <- paste0(substr(mut_snv$alteration,1,1),
                             ">",
                             substr(mut_snv$alteration,2,2))
  mut_snv$final_mut <- paste0(mut_snv$upstream, "(", mut_snv$var_type, ")", mut_snv$downstream)
  mut_snv$final_mut_1 <- paste0(substr(mut_snv$upstream, 3, 3), "(", mut_snv$var_type, ")", 
                                substr(mut_snv$downstream,1,1))
  mut_snv$final_mut_2 <- paste0(substr(mut_snv$upstream,2,2), "N(", mut_snv$var_type, ")N", 
                                substr(mut_snv$downstream,2,2))
  mut_snv$final_mut_3 <- paste0(substr(mut_snv$upstream,1,1), "NN(", mut_snv$var_type, ")NN", 
                                substr(mut_snv$downstream, 3,3))
  
  mut_snv$mut_1536 <- paste0(substr(mut_snv$upstream,2,3), "(",mut_snv$var_type,")", 
                             substr(mut_snv$downstream,1,2))
  
  return(mut_snv)
}

## Convert to maf-like structure and process mutational matrix from PCAWG
sigpromat_to_maf <- function(mat){
  maf_temp = NULL
  for(j in 3:ncol(mat)){
    temp_mat <- mat[mat[j] > 0,]
    for(i in 1:nrow(temp_mat)){
      k = temp_mat[i,j]
      maf_temp = rbind(maf_temp,
                       cbind(rep(temp_mat[i,1],k),rep(temp_mat[i,2],k),rep(colnames(temp_mat)[j],k)))
    }
  }
  return(maf_temp)
}

process_unfold_mat <- function(mat){
  mat <- as.data.frame(mat)
  colnames(mat) <- c("var_type_ori", "pentanucleotide", "cancer_sample")
  
  mat <- mat %>% 
    mutate(cancer = gsub("[::].*", "", cancer_sample)) %>%
    mutate(sample = gsub(".*::", "", cancer_sample)) %>%
    mutate(downstream = substr(pentanucleotide, 4, 5)) %>% relocate(downstream, 1) %>%
    mutate(upstream = substr(pentanucleotide, 1, 2)) %>% relocate(upstream, 1) %>%
    mutate(var_type = paste0(substr(var_type_ori, 1,1), "-", substr(var_type_ori,2,3))) %>% relocate(var_type, 1)
  return(mat)
}

## MATCH DOMINANT SIGNATURE TO SAMPLE
process_sig <- function(sig_mat, mut_mat = NULL, ngs){

  #filter by number of mutations attributed to signatures (20 in WES, 100 in WGS)
  #filter by accuracy (>=0.7)
  if(ngs == "wes"){
    sig_rel <- sig_mat[rowSums(sig_mat[,c(-1,-2,-3)])>20 & sig_mat$Accuracy > 0.7,] 
  }
  else if(ngs == "wgs"){
    sig_rel <- sig_mat[rowSums(sig_mat[,c(-1,-2,-3)])>100 & sig_mat$Accuracy > 0.7,] 
  }
  #turn to relative cont.
  for(i in 1:nrow(sig_rel)){
    sig_rel[i,c(-1,-2,-3)] <- sig_rel[i,c(-1,-2,-3)]/sum(sig_rel[i,c(-1,-2,-3)])
  }
  #select only relevant signatures
  sig_rel_F <- sig_rel %>% dplyr::select(`Cancer Types`, `Sample Names`, Accuracy,
                                         SBS1, SBS5, SBS40, SBS2, SBS13, SBS3, SBS4, 
                                         SBS6,SBS14,SBS15,SBS20,SBS21,SBS26,SBS44,
                                         SBS7a,SBS7b,SBS7c,SBS7d, SBS38, SBS10a, SBS10b, 
                                         SBS27, SBS45)
  #create features based on the signatures' etiology summed together
  sig_rel_F <- sig_rel_F %>% 
    mutate(Clock_SBS1 = SBS1) %>% 
    mutate(Clock_SBS5 = SBS5 + SBS40) %>% 
    mutate(APOBEC = SBS2+SBS13) %>% 
    mutate(HRD = SBS3) %>% 
    mutate(Tobbaco = SBS4) %>% 
    mutate(MMR = SBS6+SBS14+SBS15+SBS20+SBS21+SBS26+SBS44) %>% 
    mutate(UV = SBS7a+SBS7b+SBS7c+SBS7d+SBS38) %>% 
    mutate(POLE = SBS10a + SBS10b) %>%
    mutate(art_SBS45 = SBS45)
  
  #create target label 
  sig_rel_F <- sig_rel_F %>% mutate(sig = "other")
  for(i in 1:nrow(sig_rel_F)){
    #APOBEC
    if(sig_rel_F$APOBEC[i] > 0.3 & sig_rel_F$HRD[i]<0.4 & sig_rel_F$Tobbaco[i]<0.2 & sig_rel_F$MMR[i]<0.3){sig_rel_F$sig[i] <- "APOBEC"}
    
    else if(sig_rel_F$HRD[i] > 0.5){sig_rel_F$sig[i] <- "HRD"}
    else if(sig_rel_F$Tobbaco[i] > 0.3){sig_rel_F$sig[i] <- "Tobacco"}
    else if(sig_rel_F$MMR[i] > 0.3){sig_rel_F$sig[i] <- "MMR"}
    else if(sig_rel_F$UV[i] > 0.3){sig_rel_F$sig[i] <- "UV"}
    else if(sig_rel_F$POLE[i] > 0.2){sig_rel_F$sig[i] <- "POLE"}
    else if(sig_rel_F$Clock_SBS1[i] > 0.4){sig_rel_F$sig[i] <- "Clock_SBS1"}
    else if(sig_rel_F$Clock_SBS5[i] > 0.4){sig_rel_F$sig[i] <- "Clock_SBS5"}
  }
  if(!is.null(mut_mat)){
    mut_mat <- mut_mat %>% inner_join(sig_rel_F[,c("Sample Names","sig")], by = c("sample"="Sample Names"))
  }
  return(list(mut_mat, sig_rel_F))
}


#######----
# Efficient Generalized prediction function - VAS : 11/03/2024

cos_pred_fun_efficient <- function(sig_names, sample_names, mut, emb_len = 200, emb,
                                   creating_muts = FALSE, vep = TRUE, sim_num=0.6){

  if(creating_muts == TRUE){
    mut$mut_1536 <- paste0(mut$upstream, "(",mut$var_type,")",
                           mut$downstream)
  }

  preds <- matrix(ncol = length(sig_names),
                  nrow = length(sample_names),
                  dimnames = list(sample_names,
                                  sig_names))
  preds <- as.data.frame(cbind(sample = rownames(preds), preds))

  preds_embeddings <- mut[,colnames(mut) %in% c("sample", "Tumor_Sample_Barcode", "mut_1536")] %>%
    left_join(emb %>% rownames_to_column("mut_1536"),
              by = "mut_1536")


  if(vep == TRUE){
    preds_embeddings <- preds_embeddings %>%
      dplyr::select(!mut_1536) %>%
      group_by(Tumor_Sample_Barcode) %>%
      summarise_all(mean)

    #sort by sample name
    preds <- preds[preds_embeddings$Tumor_Sample_Barcode, ]  }

  else if(vep == FALSE){
    preds_embeddings <- preds_embeddings %>%
      dplyr::select(!mut_1536) %>%
      group_by(sample) %>%
      summarise_all(mean)

    #sort by sample name
    preds <- preds[preds_embeddings$sample, ]
  }


  for(k in 1:length(sig_names)){
    print(c("Step 1",k,"out of",length(sig_names))) ###CHANGE###
    preds[, k+1] = apply(preds_embeddings[,-1], 1,
                         FUN = MutationalPatterns::cos_sim,
                         y=as.numeric(emb[sig_names[k],]))
  }


  for(i in 2:ncol(preds)){
    print(c("Step 2",i,"out of",ncol(preds))) ###CHANGE###
    preds[,i] <- as.numeric(preds[,i])
  }

  preds_labels <- preds
  
  
  for(i in 1:nrow(preds_labels)){
    print(c("Step 3",i,"out of",nrow(preds_labels))) ###CHANGE###
    preds_labels[i, 2:ncol(preds_labels)] <-
      ifelse(as.numeric(preds[i,2:ncol(preds)]) >= 
               max(as.numeric(preds[i,2:ncol(preds)])) &
               max(as.numeric(preds[i,2:ncol(preds)])) > sim_num,
             "Yes","No")
  }
  
  return(list(preds, preds_embeddings, preds_labels))
}


###final pred function
confmat_preds <- function(preds, sig_col=11, rounds=FALSE, filt = FALSE, cutoff = 2){
  if (rounds == TRUE){
    for(i in 1:nrow(preds)){
      preds[i,2:10] <- ifelse(preds[i,2:10] == max(as.numeric(preds[i,2:10])) & max(as.numeric(preds[i,2:10])) > cutoff,"Yes","No")
    }
  }
  preds$pred <- "No prediction"
  for(i in 1:nrow(preds)){
    for(j in 2:10){
      if(preds[i,j] == "Yes"){
        preds$pred[i] <- colnames(preds)[j]
      }
    }
  }
  
  if(filt ==TRUE){
    preds <- preds %>% filter(!sig.y %in% c("other","HRD")) %>% filter(!pred %in% c("other","HRD"))
  }
  
  return(preds)
}

#for unlabeled data
add_pred <- function(preds){
  preds$pred <- "No prediction"
  for(i in 1:nrow(preds)){
    for(j in 2:10){
      if(preds[i,j] == "Yes"){
        preds$pred[i] <- colnames(preds)[j]
      }
    }
  }
  return(preds)
}

## add cancer type MSK based panels
can_type_add <- function(preds, clin){
  
  preds <- preds %>% left_join(clin[,c("CANCER_TYPE","SAMPLE_ID")], by=c("sample"="SAMPLE_ID"))
  preds_per_cancer <- as.data.frame(table(preds$CANCER_TYPE, preds$pred))
  colnames(preds_per_cancer) <- c("Cancer Type", "Sig", "Freq")
  
  return(list(preds,preds_per_cancer))
}


## num of type sig appear in cancer type - msk
num_per_cancer <- function(msk_preds, cancer_type, sig){
  s <- sum(msk_preds[msk_preds$CANCER_TYPE == cancer_type, sig] == "Yes")
  return(s)
}


####----
# benchmark nnls and mix functions
preds_cancer_type_nnls_onco_noMutCat <- function(preds_joined, sigs, num_above){
  samples <- NULL
  
  if(length(sigs) > 1){
    preds_joined$sigs_cont <- 0
    for(x in sigs){
      preds_joined$sigs_cont <- preds_joined$sigs_cont + preds_joined[,x]
    }
    
    samples <- unique(c(samples, preds_joined$sample[preds_joined$sigs_cont > num_above]))
  }
  
  else{
    for(x in sigs){
      samples <- unique(c(samples, preds_joined$sample[preds_joined[,x] > num_above]))
    }
  }
  return(table(preds_joined$ONCOTREE[preds_joined$sample %in%
                                       samples]))
}


preds_cancer_type_nnls_onco <- function(preds_joined, sigs, num_above, mutcat){
  samples <- NULL
  
  preds_joined <- preds_joined %>%
    filter(mutCat == mutcat)
  
  if(length(sigs) > 1){
    preds_joined$sigs_cont <- 0
    for(x in sigs){
      preds_joined$sigs_cont <- preds_joined$sigs_cont + preds_joined[,x]
    }
    
    samples <- unique(c(samples, preds_joined$sample[preds_joined$sigs_cont > num_above]))
  }
  
  else{
    for(x in sigs){
      samples <- unique(c(samples, preds_joined$sample[preds_joined[,x] > num_above]))
    }
  }
  return(table(preds_joined$ONCOTREE[preds_joined$sample %in%
                                       samples]))
}



###############################
#       50% rel for benchmark
##############################
process_sig_benchm <- function(processed_sig_mat, siglabels){
  temp <- NULL
  for(x in siglabels){
    if(x != "Clock_SBS5"){
      temp <- rbind(
        temp,
        processed_sig_mat[processed_sig_mat[,x] >= 0.5,])
    }
    else if(x == "Clock_SBS5"){
      temp <- rbind(
        temp,
        processed_sig_mat[processed_sig_mat[,x] >= 0.75,])
    }

  }
  return(temp)
}


#######################
# create mean + se + ci for sampling experiments
######################

mean_per_experiment <- function(df){
  df_mean <- df %>%
    group_by(sig, mut_num) %>%
    dplyr::summarise(
      mean_sens = mean(Sensitivity),
      sd_sens = sd(Sensitivity),
      n_sens = n(),
      se_sens = sd_sens / sqrt(n_sens),
      
      mean_ppv = mean(`Pos Pred Value`),
      sd_ppv = sd(`Pos Pred Value`),
      n_ppv = n(),
      se_ppv = sd_ppv / sqrt(n_ppv),       
      
      mean_spec = mean(Specificity),
      sd_spec = sd(Specificity),
      n_spec = n(),
      se_spec = sd_spec / sqrt(n_spec),
      
      mean_npv = mean(`Neg Pred Value`),
      sd_npv = sd(`Neg Pred Value`),
      n_npv = n(),
      se_npv = sd_npv / sqrt(n_npv)     
    )
  
  df_mean$ci_sens <- df_mean$se_sens * qt(.975, df_mean$n_sens - 1)
  df_mean$ci_ppv <- df_mean$se_ppv * qt(.975, df_mean$n_ppv - 1)
  df_mean$ci_spec <- df_mean$se_spec * qt(.975, df_mean$n_spec - 1)
  df_mean$ci_npv <- df_mean$se_npv * qt(.975, df_mean$n_npv - 1)

  
  df_mean$sig <- gsub("Class: ", "", df_mean$sig)
  
  return(df_mean)
}

