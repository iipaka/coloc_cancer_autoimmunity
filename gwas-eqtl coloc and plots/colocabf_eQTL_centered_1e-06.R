# ==============================================================================
# Purpose: Master loop for coloc analysis across all target regions
# ==============================================================================

library(coloc)
library(dplyr)

p_threshold <- "1e-06"

sig_snps <- read.csv("Output/recentered_targets.csv")
sig_snps$clean_snp <- sub(":.*", "", sig_snps$New_Centre_SNP)

gwas_info <- read.csv("Data_raw/gwas_summarystats.csv")

eqtl_files <- list.files("Data_processed", pattern = "^eqtl_.*_centered\\.rds$", full.names = TRUE)
message(paste("\nFound", length(eqtl_files), "total centered eQTL files to test."))

results_list <- list()

# ==============================================================================
# LOOP THROUGH SNPs
# ==============================================================================

results_counter <- 1 

for(i in seq_along(eqtl_files)) {
  
  eqtl_f <- eqtl_files[i]
  clean_name <- gsub("_centered\\.rds$", "", basename(eqtl_f))
  parts <- strsplit(clean_name, "_")[[1]]
  
  rsid <- parts[2]
  gene_name <- parts[3]
  
  target_rows <- sig_snps[sig_snps$New_Centre_SNP == rsid, ]
  
  unique_traits <- unique(target_rows$Stronger_Trait)
  
  for(trait_name in unique_traits) {
    
    # Isolate the coordinate data for this specific unique trait
    trait_info <- target_rows[target_rows$Stronger_Trait == trait_name, ][1, ]
    target_chr <- trait_info$CHR
    target_pos <- trait_info$New_Centre_POS
    
    message(sprintf("[%d/%d] Testing: %s vs %s (GWAS: %s)", i, length(eqtl_files), rsid, gene_name, trait_name))
    
    trait_gcst <- gwas_info$name[gwas_info$Disease == trait_name]
    file_gwas <- file.path("Data_processed", paste0(rsid, "_", trait_gcst, "_", p_threshold, "_centered.rds"))
    
    list_eqtl <- readRDS(eqtl_f)
    list_gwas <- readRDS(file_gwas)
    
    common_snps <- intersect(list_gwas$snp, list_eqtl$snp)
    if(length(common_snps) == 0) {
      message(sprintf("  -> WARNING: 0 overlapping SNPs for %s vs %s. Skipping...", gene_name, trait_name))
      next
    }
    
    idx_gwas <- which(list_gwas$rsid == rsid)
    idx_eqtl <- which(list_eqtl$rsid == rsid)
    beta_gwas <- list_gwas$beta[idx_gwas[1]]
    beta_eqtl <- list_eqtl$beta[idx_eqtl[1]]
    
    coloc_res <- coloc.abf(dataset1 = list_gwas, dataset2 = list_eqtl)
    results <- coloc_res$summary
    
    results_list[[results_counter]] <- data.frame(
      Target_SNP = rsid,
      CHR = target_chr,
      POS = target_pos,
      GWAS_Trait = trait_name,
      eQTL_Gene = gene_name,
      GWAS_Beta = beta_gwas, 
      eQTL_Beta = beta_eqtl,
      nsnps = results["nsnps"],
      PP_H0 = results["PP.H0.abf"],
      PP_H1 = results["PP.H1.abf"],
      PP_H2 = results["PP.H2.abf"],
      PP_H3 = results["PP.H3.abf"],
      PP_H4 = results["PP.H4.abf"]
    )
    
    results_counter <- results_counter + 1
  }
}

# ==============================================================================
# COMPILE AND SAVE RESULTS
# ==============================================================================

final_results <- bind_rows(results_list) %>%
  arrange(Target_SNP, desc(PP_H4))

output_file <- file.path("Output", paste0("coloc_results_eqtl_gwas_summary_", p_threshold, "_centered.csv"))
write.csv(final_results, output_file, row.names = FALSE)