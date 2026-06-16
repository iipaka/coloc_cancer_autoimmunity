# ==============================================================================
# Purpose: Extract shifted 1MB eQTLGen windows centered on the true regional peak
# ==============================================================================

library(glue)
library(data.table)
library(dplyr)

options(warn = 1)
window <- 1000000 / 2 

targets <- fread("Output/recentered_targets.csv")
eqtl_file <- "/rds/user/xz418/hpc-work/Data_eQTLGen/Blood_eQTLGen_hg38_sorted.tsv.gz"

tmp_dir <- "tmp_dir/"
output_dir <- "Data_processed/"

eqtl_cols <- c("CHR38", "BP38", "SNPID", "CHR19", "BP19", "REF", "ALT", "ALT_FREQ", 
               "P", "Zscore", "Gene_id", "Gene_name", "GeneChr", "GenePos", 
               "NrCohorts", "NrSamples", "FDR", "BonferroniP", "BETA", "SE")

for(k in 1:nrow(targets)) {
  
  t_snp <- targets$New_Centre_SNP[k] 
  t_chr <- targets$CHR[k] 
  t_pos <- targets$New_Centre_POS[k] 
  
  pos_min <- t_pos - window
  pos_max <- t_pos + window
  
  tabix_cmd <- glue("tabix {eqtl_file} {t_chr}:{pos_min}-{pos_max}")
  
  region_data <- tryCatch(
    fread(cmd = tabix_cmd, col.names = eqtl_cols, tmpdir = tmp_dir),
    error = function(e) return(NULL)
  )
  
  region_data <- region_data[!is.na(BETA) & !is.na(SE) & !is.na(P)]
  region_data[, snp_id := paste0(CHR38, ":", BP38)]
  region_data <- region_data[!duplicated(paste0(snp_id, "_", Gene_name))]
  
  window_genes <- unique(region_data$Gene_name)
  message(paste0(" -> Found ", length(window_genes), " unique genes. Generating files..."))
  
  for(current_gene in window_genes) {
    gene_data <- region_data[Gene_name == current_gene]
    
    coloc_input <- list(
      snp = gene_data$snp_id, 
      rsid = gene_data$SNPID, 
      chr = t_chr,
      beta = gene_data$BETA, 
      varbeta = gene_data$SE^2, 
      pvalues = gene_data$P,
      position = gene_data$BP38, 
      ref = gene_data$REF, 
      alt = gene_data$ALT,
      type = "quant", 
      N = gene_data$NrSamples,
      MAF = ifelse(gene_data$ALT_FREQ > 0.5, 1 - gene_data$ALT_FREQ, gene_data$ALT_FREQ) 
    )
    
    save_name <- paste0(output_dir, "eqtl_", t_snp, "_", current_gene, "_centered.rds")
    saveRDS(coloc_input, save_name)
  }
}