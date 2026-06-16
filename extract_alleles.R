# ==============================================================================
# Purpose: Check ref/alt allele mismatch
# ==============================================================================

library(dplyr)

p_threshold <- "1e-06"
out_file <- "Output/Allele_refalt_check.csv"

gwas_info <- read.csv("Data_raw/gwas_summarystats.csv")
gwas_info$Disease <- trimws(gwas_info$Disease)

targets <- read.csv(paste0("Output/sig_snps_threshold_", p_threshold, ".csv"))

results_list <- list()

for(i in 1:nrow(targets)) {
  
  rsid <- targets$rsID[i] 
  trait1_name <- trimws(targets$Cancer_Type[i])
  trait2_name <- trimws(targets$Name[i])
  
  message(sprintf("[%d/%d] Extracting %s (%s vs %s)", i, nrow(targets), rsid, trait1_name, trait2_name))
  
  trait1_gcst <- gwas_info$name[gwas_info$Disease == trait1_name]
  trait2_gcst <- gwas_info$name[gwas_info$Disease == trait2_name]
  
  file_trait1 <- file.path("Data_processed", paste0(rsid, "_", trait1_gcst, "_", p_threshold, ".rds"))
  file_trait2 <- file.path("Data_processed", paste0(rsid, "_", trait2_gcst, "_", p_threshold, ".rds"))
  
  d_t1 <- readRDS(file_trait1)
  d_t2 <- readRDS(file_trait2)
  df1 <- as.data.frame(d_t1)
  df2 <- as.data.frame(d_t2)
  
  row1 <- df1[grepl(rsid, df1$rsid), c("snp", "rsid", "ref", "alt", "beta")]
  row2 <- df2[grepl(rsid, df2$rsid), c("snp", "rsid", "ref", "alt", "beta")]
  
  if(nrow(row1) > 0) {
    row1 <- row1[1, ] 
  } else {
    message("  -> WARNING: SNP missing in Trait 1. Imputing NAs...")
    row1 <- data.frame(snp=NA, rsid=rsid, ref=NA, alt=NA, beta=NA)
  }
  
  if(nrow(row2) > 0) {
    row2 <- row2[1, ] 
  } else {
    message("  -> WARNING: SNP missing in Trait 2. Imputing NAs...")
    row2 <- data.frame(snp=NA, rsid=rsid, ref=NA, alt=NA, beta=NA)
  }
  
  combined_row <- data.frame(
    Target_SNP = rsid,
    Trait1 = trait1_name,
    Trait1_Ref = row1$ref,
    Trait1_Alt = row1$alt,
    Trait1_Beta = row1$beta,
    Trait2 = trait2_name,
    Trait2_Ref = row2$ref,
    Trait2_Alt = row2$alt,
    Trait2_Beta = row2$beta
  )
  
  results_list[[i]] <- combined_row
}

final_results <- bind_rows(results_list)
write.csv(final_results, out_file, row.names = FALSE)

message("\nExtraction complete! Results saved to: ", out_file)