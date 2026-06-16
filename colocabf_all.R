# ==============================================================================
# Purpose: Master loop for coloc analysis across all target regions
# ==============================================================================

library(coloc)
library(dplyr)

# 1. Load reference files 
sig_snps <- read.csv("Output/sig_snps_threshold_5e-08.csv")
gwas_info <- read.csv("Data_raw/gwas_summarystats.csv")

# Initialize an empty list to store results
results_list <- list()

# ==============================================================================
# 2. LOOP THROUGH SNPs
# ==============================================================================

for(i in 1:nrow(sig_snps)) {
  
  # A. Extract relevant columns for mapping
  rsid        <- sig_snps[i, 1]
  target_chr  <- sig_snps[i, 2]
  target_pos  <- sig_snps[i, 3]
  cancer_name <- sig_snps[i, 16]
  immune_name <- sig_snps[i, 10]
  
  message(paste0("\n[", i, "/", nrow(sig_snps), "] Running coloc for ", rsid, " (", cancer_name, " vs ", immune_name, ")"))
  
  # B. Dynamically find the correct GCST IDs 
  cancer_gcst <- gwas_info$name[gwas_info$Disease == cancer_name]
  immune_gcst <- gwas_info$name[gwas_info$Disease == immune_name]
  
  # C. Construct the exact file paths
  file_cancer <- paste0("Data_processed/", rsid, "_", cancer_gcst, ".rds")
  file_immune <- paste0("Data_processed/", rsid, "_", immune_gcst, ".rds")
  
  # D. Load the processed lists
  list_cancer <- readRDS(file_cancer)
  list_immune <- readRDS(file_immune)
  
  # E. Run coloc
  coloc_res <- coloc.abf(dataset1 = list_cancer, dataset2 = list_immune)
  
  # F. Results
  results <- coloc_res$summary
  
  # I. Store the results in a dataframe row
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
# 3. COMPILE AND SAVE RESULTS
# ==============================================================================

final_results <- bind_rows(results_list)

output_file <- "Output/coloc_results_summary.csv"
write.csv(final_results, output_file, row.names = FALSE)

message(paste0("\nAll coloc analyses complete. Summary saved to: ", output_file))