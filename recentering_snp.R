# ==============================================================================
# Purpose: Find the regional lead SNP in the Stronger Trait window
# ==============================================================================

library(dplyr)

p_threshold <- "1e-06"
targets_file <- paste0("Output/eQTL_coloc_gwas_targets_", p_threshold, ".csv")
out_file <- "Output/recentered_targets.csv"

targets <- read.csv(targets_file)
gwas_info <- read.csv("Data_raw/gwas_summarystats.csv")
gwas_info$Disease <- gwas_info$Disease

targets$New_Centre_SNP <- NA
targets$New_Centre_POS <- NA
targets$New_Centre_P <- NA

for(i in 1:nrow(targets)) {
  
  old_snp <- targets$SNP[i]
  s_trait <- targets$Stronger_Trait[i]
  
  gcst <- gwas_info$name[gwas_info$Disease == s_trait]
  
  rds_file <- file.path("Data_processed", paste0(old_snp, "_", gcst, "_", p_threshold, ".rds"))

  d <- readRDS(rds_file)
  df <- as.data.frame(d)
    
  best_row <- df[which.min(df$pvalues), ]
    
  targets$New_Centre_SNP[i] <- best_row$rsid
  targets$New_Centre_POS[i] <- best_row$position
  targets$New_Centre_P[i] <- best_row$pvalues
    
}

write.csv(targets, out_file, row.names = FALSE)