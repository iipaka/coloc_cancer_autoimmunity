# ==============================================================================
# Purpose: Extract shifted 1MB GWAS windows centered on the true regional peak
# ==============================================================================

library(glue)
library(data.table)
library(dplyr)

options(warn = 1)
p_threshold <- "1e-06"
window <- 1000000 / 2

targets <- fread("Output/recentered_targets.csv")
gwas_summarystats <- fread("Data_raw/gwas_summarystats.csv")
gwas_summarystats$Disease <- gwas_summarystats$Disease

tmp_dir <- "tmp_dir/"
output_dir <- "Data_processed/"
dir.create(tmp_dir, showWarnings = FALSE)

for(k in 1:nrow(targets)) {
  
  t_snp <- targets$New_Centre_SNP[k]
  t_chr <- targets$CHR[k]
  t_pos <- targets$New_Centre_POS[k]
  s_trait <- targets$Stronger_Trait[k]
  
  gwas_row <- gwas_summarystats[Disease == s_trait]
  
  f_path <- gwas_row$file_path[1]
  s_name <- gwas_row$name[1]
  s_N <- gwas_row$N[1]
  s_s <- gwas_row$s[1]
  
  pos_min <- t_pos - window
  pos_max <- t_pos + window
  
  header_cmd <- paste("zcat", f_path, "| head -n 1") 
  header <- names(fread(cmd = header_cmd)) 
  
  if("hm_chrom" %in% header) {
    col_idx_chr <- 3; col_idx_pos <- 4
    old_cols <- c("hm_rsid", "hm_chrom", "hm_pos", "hm_odds_ratio", "p_value", "hm_beta", "standard_error", "hm_effect_allele", "hm_other_allele")
    new_cols <- c("rsid", "CHR", "POS", "OR", "P", "BETA", "SE", "ALT", "REF")
  } else if ("CHR38" %in% header) {
    col_idx_chr <- which(header == "CHR38"); col_idx_pos <- which(header == "BP38")
    old_cols <- c("SNPID", "CHR38", "BP38", "P", "BETA", "SE", "ALT", "REF")
    new_cols <- c("rsid", "CHR", "POS", "P", "BETA", "SE", "ALT", "REF")
  } else if ("CHR" %in% header) {
    col_idx_chr <- which(header == "CHR"); col_idx_pos <- which(header == "POS")
    old_cols <- c("snp", "CHR", "POS", "P", "BETA", "SE", "ALT", "REF")
    new_cols <- c("rsid", "CHR", "POS", "P", "BETA", "SE", "ALT", "REF")
  } else {
    col_idx_chr <- 1; col_idx_pos <- 2
    old_cols <- c("rsid", "chromosome", "base_pair_location", "p_value", "beta", "standard_error", "effect_allele", "other_allele")
    new_cols <- c("rsid", "CHR", "POS", "P", "BETA", "SE", "ALT", "REF")
  }
  
  regionfetch_cmd <- glue("zcat {f_path} | awk 'NR==1 || (${col_idx_chr}=={t_chr} && ${col_idx_pos} >= {pos_min} && ${col_idx_pos} <= {pos_max})'")
  region_data <- fread(cmd=regionfetch_cmd, tmpdir=tmp_dir)
  
  setnames(region_data, old=old_cols, new=new_cols, skip_absent = TRUE)
  if("BETA" %in% names(region_data)) region_data[, BETA := as.numeric(BETA)]
  if("SE" %in% names(region_data)) region_data[, SE := as.numeric(SE)]
  if("P" %in% names(region_data)) region_data[, P := as.numeric(P)]
  
  if ("BETA" %in% names(region_data) && "OR" %in% names(region_data)) {
    region_data[is.na(BETA) & !is.na(OR), BETA := log(as.numeric(OR))]
  } else if (!"BETA" %in% names(region_data) && "OR" %in% names(region_data)) {
    region_data[, BETA := log(as.numeric(OR))]
  }
  
  region_data <- region_data[!is.na(BETA) & !is.na(SE)]
  if(min(region_data$P, na.rm=TRUE) == 0) region_data[P==0, P := exp(pnorm(-abs(BETA)/SE, log.p=TRUE))]
  
  valid_alleles <- c("A", "C", "G", "T")
  region_data <- region_data[REF %in% valid_alleles & ALT %in% valid_alleles]
  region_data[, snp := paste0(CHR, ":", POS)]
  region_data <- region_data[!duplicated(snp)]
  
  coloc_input <- list(
    snp=region_data$snp, 
    rsid=region_data$rsid, 
    chr=t_chr, 
    beta=region_data$BETA, 
    varbeta=region_data$SE^2, 
    pvalues=region_data$P, 
    position=region_data$POS, 
    ref=region_data$REF, 
    alt=region_data$ALT, 
    type="cc", 
    s=s_s, 
    N=s_N
  )
  
  save_name <- paste0(output_dir, t_snp, "_", s_name, "_", p_threshold, "_centered.rds")
  saveRDS(coloc_input, save_name)
}
