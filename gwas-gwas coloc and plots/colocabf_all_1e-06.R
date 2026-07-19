# ==============================================================================
# Purpose: Master loop for coloc analysis across all target regions
# ==============================================================================

library(coloc)
library(dplyr)

p_threshold <- 1e-6

sig_snps <- read.csv(file.path("Output", paste0("sig_snps_threshold_", p_threshold, ".csv")))
gwas_info <- read.csv(file.path("Data_raw", "gwas_summarystats.csv"))

results_list <- list()

# ==============================================================================
# LOOP THROUGH SNPs
# ==============================================================================

for(i in 1:nrow(sig_snps)) {
  
  #Extract relevant columns for mapping
  rsid <- sig_snps[i, 1]
  target_chr <- sig_snps[i, 2]
  target_pos <- sig_snps[i, 3]
  cancer_name <- sig_snps[i, 16]
  immune_name <- sig_snps[i, 10]
  
  message(paste0("\n[", i, "/", nrow(sig_snps), "] Running coloc for ", rsid, " (", cancer_name, " vs ", immune_name, ")"))
  
  #Find the correct GCST IDs 
  cancer_gcst <- gwas_info$name[gwas_info$Disease == cancer_name]
  immune_gcst <- gwas_info$name[gwas_info$Disease == immune_name]
  
  #Construct file paths
  file_cancer <- file.path("Data_processed", paste0(rsid, "_", cancer_gcst, "_", p_threshold, ".rds"))
  file_immune <- file.path("Data_processed", paste0(rsid, "_", immune_gcst, "_", p_threshold, ".rds"))
  
  #Load rds
  list_cancer <- readRDS(file_cancer)
  list_immune <- readRDS(file_immune)
  
  #Run coloc
  coloc_res <- coloc.abf(dataset1 = list_cancer, dataset2 = list_immune)
  results <- coloc_res$summary
  
  #Store results in a dataframe row
  res_row <- data.frame(
    SNP = rsid,
    CHR = target_chr,
    POS = target_pos,
    Cancer_trait = cancer_name,
    Immune_trait = immune_name,
    nsnps = results["nsnps"],
    PP_H0 = results["PP.H0.abf"],
    PP_H1 = results["PP.H1.abf"],
    PP_H2 = results["PP.H2.abf"],
    PP_H3 = results["PP.H3.abf"],
    PP_H4 = results["PP.H4.abf"]
  )
  
  results_list[[i]] <- res_row
  message(sprintf("  -> Done. PP.H4: %.4f | PP.H3: %.4f", results["PP.H4.abf"], results["PP.H3.abf"]))
}

# ==============================================================================
# COMPILE AND SAVE RESULTS
# ==============================================================================

final_results <- bind_rows(results_list)

output_file <- file.path("Output", paste0("coloc_results_summary_", p_threshold, ".csv"))
write.csv(final_results, output_file, row.names = FALSE)

message(paste0("\nAll coloc analyses complete. Summary saved to: ", output_file))