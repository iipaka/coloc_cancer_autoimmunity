# ==============================================================================
# Purpose: Master loop for coloc analysis across all target regions
# ==============================================================================

library(coloc)
library(dplyr)

# 1. Load reference files 
p_threshold <- 1e-6

sig_snps <- read.csv(file.path("Output", paste0("eQTL_coloc_gwas_targets_", p_threshold, ".csv")))
gwas_info <- read.csv(file.path("Data_raw", "gwas_summarystats.csv"))

eqtl_files <- list.files("Data_processed", pattern = "^eqtl_.*\\.rds$", full.names = TRUE)
message(paste("\nFound", length(eqtl_files), "total eQTL files to test."))

results_list <- list()

# ==============================================================================
# LOOP THROUGH SNPs
# ==============================================================================

for(i in seq_along(eqtl_files)) {
  
  eqtl_f <- eqtl_files[i]
  
  clean_name <- gsub("\\.rds$", "", basename(eqtl_f))
  parts <- strsplit(clean_name, "_")[[1]]
  
  rsid <- parts[2]
  gene_name <- parts[3]
  
  target_row <- sig_snps[sig_snps[, 1] == rsid, ][1, ]
  target_chr <- target_row[, 2]
  target_pos <- target_row[, 3]
  trait_name <- target_row[, 10]
  
  message(sprintf("[%d/%d] Testing: %s vs %s (GWAS: %s)", i, length(eqtl_files), rsid, gene_name, trait_name))
  
  trait_gcst <- gwas_info$name[gwas_info$Disease == trait_name]
  file_gwas <- file.path("Data_processed", paste0(rsid, "_", trait_gcst, "_", p_threshold, ".rds"))
  
  
  list_eqtl <- readRDS(eqtl_f)
  list_gwas <- readRDS(file_gwas)
  
  common_snps <- intersect(list_gwas$snp, list_eqtl$snp)
  if(length(common_snps) == 0) {
    message(sprintf("  -> WARNING: 0 overlapping SNPs for %s. Skipping...", gene_name))
    next
  }
  
  coloc_res <- coloc.abf(dataset1 = list_gwas, dataset2 = list_eqtl)

  results <- coloc_res$summary
  
  results_list[[i]] <- data.frame(
    Target_SNP = rsid,
    CHR = target_chr,
    POS = target_pos,
    GWAS_Trait = trait_name,
    eQTL_Gene = gene_name,
    nsnps = results["nsnps"],
    PP_H0 = results["PP.H0.abf"],
    PP_H1 = results["PP.H1.abf"],
    PP_H2 = results["PP.H2.abf"],
    PP_H3 = results["PP.H3.abf"],
    PP_H4 = results["PP.H4.abf"]
  )
}

# ==============================================================================
# COMPILE AND SAVE RESULTS
# ==============================================================================

final_results <- bind_rows(results_list) %>%
  arrange(Target_SNP, desc(PP_H4))

output_file <- file.path("Output", paste0("coloc_results_eqtl_gwas_summary_", p_threshold, ".csv"))
write.csv(final_results, output_file, row.names = FALSE)

message(paste0("\nSuccess! Summary saved to: ", output_file))