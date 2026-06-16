# ==============================================================================
# Purpose: Extract and clean specific genomic regions for coloc analysis
# ==============================================================================

library(glue)
library(data.table)
library(dplyr)
library(R.utils)

# ------------------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------------------
options(warn = 1)

p_threshold<-1e-6

candidate_snps <- fread(paste0("Output/sig_snps_threshold_", p_threshold, ".csv"))
gwas_summarystats <- fread("Data_raw/gwas_summarystats.csv")

tmp_dir <- "tmp_dir/"
output_dir <- "Data_processed/"

window <- 1000000/2

for(i in 1:nrow(gwas_summarystats)) {
    s_name <- gwas_summarystats$name[i]
    f_path <- gwas_summarystats$file_path[i]
    s_N <- gwas_summarystats$N[i]
    s_s <- gwas_summarystats$s[i]
    
    message(paste0("\nProcessing Study", s_name))
    
    header_cmd <- paste("zcat", f_path, "| head -n 1")
    header <- names(fread(cmd = header_cmd))
    
    if("hm_chrom" %in% header) {
      col_idx_chr <- 3
      col_idx_pos <- 4
      
      old_cols <- c("hm_rsid", "hm_chrom", "hm_pos", "hm_odds_ratio", "p_value", 
                    "hm_beta", "standard_error", "hm_effect_allele", "hm_other_allele",
                    "hm_ci_lower", "hm_ci_upper")
      new_cols <- c("rsid", "CHR", "POS", "OR", "P", "BETA", "SE", "ALT", "REF", "CI_LOWER", "CI_UPPER")
      
    } else if ("CHR38" %in% header) {
      col_idx_chr <- which(header == "CHR38")
      col_idx_pos <- which(header == "BP38")
      
      old_cols <- c("SNPID", "CHR38", "BP38", "P", "BETA", "SE", "ALT", "REF")
      new_cols <- c("rsid", "CHR", "POS", "P", "BETA", "SE", "ALT", "REF")
      
    } else if ("CHR" %in% header) {
      col_idx_chr <- which(header == "CHR")
      col_idx_pos <- which(header == "POS")
      
      old_cols <- c("snp", "CHR", "POS", "P", "BETA", "SE", "ALT", "REF")
      new_cols <- c("rsid", "CHR", "POS", "P", "BETA", "SE", "ALT", "REF")
      
    } else {
      col_idx_chr <- 1
      col_idx_pos <- 2
      
      old_cols <- c("rsid", "chromosome", "base_pair_location", "p_value", 
                    "beta", "standard_error", "effect_allele", "other_allele")
      new_cols <- c("rsid", "CHR", "POS", "P", "BETA", "SE", "ALT", "REF")
    }
  
  for(k in 1:nrow(candidate_snps)){
    t_snp <- candidate_snps[[1]][k]
    t_chr <- candidate_snps[[2]][k]
    
    pos_min <- candidate_snps[[3]][k] - window
    pos_max <- candidate_snps[[3]][k] + window

# ------------------------------------------------------------------------------
# 2. SNP slicing 
# ------------------------------------------------------------------------------

    regionfetch_cmd <- glue("zcat {f_path} | awk 'NR==1 || (${col_idx_chr}=={t_chr} && ${col_idx_pos} >= {pos_min} && ${col_idx_pos} <= {pos_max})'")
    region_data <- fread(cmd=regionfetch_cmd, tmpdir=tmp_dir)

# ------------------------------------------------------------------------------
# 3. Data cleaning 
# ------------------------------------------------------------------------------

    setnames(region_data, old=old_cols, new=new_cols, skip_absent = TRUE)
  
    if("BETA" %in% names(region_data)) region_data[, BETA := as.numeric(BETA)]
    if("SE" %in% names(region_data)) region_data[, SE := as.numeric(SE)]
    if("P" %in% names(region_data)) region_data[, P := as.numeric(P)]
    
    if ("BETA" %in% names(region_data) && "OR" %in% names(region_data)) {
      region_data[is.na(BETA) & !is.na(OR), BETA := log(as.numeric(OR))]
    } else if (!"BETA" %in% names(region_data) && "OR" %in% names(region_data)) {
      region_data[, BETA := log(as.numeric(OR))]
    }
    
    if ("CI_UPPER" %in% names(region_data) && "CI_LOWER" %in% names(region_data)) {
      if ("SE" %in% names(region_data)) {
        region_data[is.na(SE) & !is.na(CI_UPPER) & !is.na(CI_LOWER), 
                    SE := (log(as.numeric(CI_UPPER)) - log(as.numeric(CI_LOWER))) / 3.92]
      } else {
        region_data[, SE := (log(as.numeric(CI_UPPER)) - log(as.numeric(CI_LOWER))) / 3.92]
      }
    }
    
    n_start <- nrow(region_data)
    n_curr  <- n_start
    message(paste0("Initial variants: ", n_start))
    
    
    region_data <- region_data[!is.na(BETA) & !is.na(SE)]
    n_new <- nrow(region_data)
    message(paste0("- Lost to Missing Stats: ", n_curr - n_new))
    n_curr <- n_new
    
    if(min(region_data$P, na.rm=TRUE) == 0) {
      region_data[P==0, P := exp(pnorm(-abs(BETA)/SE, log.p=TRUE))]
    }
    
    if(is.character(region_data$CHR)) region_data[, CHR := as.numeric(CHR)]
    
    region_data <- region_data[!is.na(CHR)]
    n_new <- nrow(region_data)
    message(paste0("- Lost to Bad CHR: ", n_curr - n_new))
    n_curr <- n_new
    
    valid_alleles <- c("A", "C", "G", "T")
    region_data <- region_data[REF %in% valid_alleles & ALT %in% valid_alleles]
    n_new <- nrow(region_data)
    message(paste0("- Lost to Non-SNVs: ", n_curr - n_new))
    n_curr <- n_new
    
    region_data <- region_data[!(BETA == Inf | BETA == -Inf | SE == Inf | SE <= 0)]
    n_new <- nrow(region_data)
    message(paste0("- Lost to Infinite/Bad SE: ", n_curr - n_new))
    n_curr <- n_new
    
    region_data[, snp := paste0(CHR, ":", POS)]
    region_data <- region_data[!duplicated(snp)]
    n_new <- nrow(region_data)
    message(paste0("- Lost to Duplicates: ", n_curr - n_new))
    n_curr <- n_new
    
    n_total_lost <- n_start - n_curr
    pct_lost <- (n_total_lost / n_start) * 100
    
    message(paste0("-----------------------------"))
    message(paste0("FINAL CLEAN COUNT: ", n_curr))
    message(paste0("TOTAL REMOVED: ", n_total_lost, " (", round(pct_lost, 2), "%)"))
    
    # --------------------------------------------------------------------------
    # p-value QC
    # --------------------------------------------------------------------------
    
    z_score_check <- region_data$BETA / region_data$SE
    p_theoretical <- 2 * pnorm(-abs(z_score_check))
    
    cor_p <- cor(-log10(region_data$P + 1e-300), 
                 -log10(p_theoretical + 1e-300), 
                 use = "complete.obs")
    
    message(paste0("P-value correlation: ", round(cor_p, 4)))
    
    if(!is.na(cor_p) && cor_p < 0.90) {
      message("Warning: P-values do not match Beta/SE. (Corr < 0.90)")
    }
    
# ------------------------------------------------------------------------------
# 4. coloc object creation
# ------------------------------------------------------------------------------

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
    
    save_name <- paste0(output_dir, t_snp, "_", s_name, "_", p_threshold, ".rds")
    saveRDS(coloc_input, save_name)
  }
}

message("Extraction complete.")