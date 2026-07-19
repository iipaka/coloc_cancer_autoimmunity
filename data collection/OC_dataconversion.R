# ==============================================================================
# Purpose: Standardise Ovarian Cancer & HGS OC GWAS
# ==============================================================================

library(data.table)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

# ------------------------------------------------------------------------------
# SETUP
# ------------------------------------------------------------------------------

test_file <- "Data_raw/Phelan_Archive/Summary_chr6.txt"

header <- names(fread(test_file, nrows=0))
print(header)

#Chr at 3rd col, pos at 4th col

dt <- fread("Data_Raw/Phelan_Archive/Summary_chr9.txt", nrows = 10)
print(dt[1:10, .(OrigSNPname, `1000G_SNPname`, Chromosome, Position)])

input_folder <- "Data_Raw/Phelan_Archive/" 

out_overall <- "Data_Raw/Ovarian_Cancer_Overall.tsv.gz"
out_hgs <- "Data_Raw/Ovarian_Cancer_HGS.tsv.gz"


# ------------------------------------------------------------------------------
# SAMPLING FOR BUILD VERIFICATION
# ------------------------------------------------------------------------------

sampling_results <- list()
files <- list.files(input_folder, pattern = "Summary_chr", full.names = TRUE)

for(f in files) {
  fname <- basename(f)
  
  temp_dt <- fread(f, select = c("OrigSNPname", "Chromosome", "Position"))
  rs_only <- temp_dt[grepl("^rs", OrigSNPname)]
  
  if(nrow(rs_only) >= 10) {
    #Randomly sample 10 rows
    sample_sub <- rs_only[sample(.N, 10)]
    sample_sub[, Source_File := fname]
    sampling_results[[fname]] <- sample_sub
  } else {
    message(paste("Warning: Fewer than 10 rsIDs found in", fname))
    rs_only[, Source_File := fname]
    sampling_results[[fname]] <- rs_only
  }
}

verification_set <- rbindlist(sampling_results)

fwrite(verification_set, "Output/SNP_build_check_sample.csv")

#Checking SNPs obtained using the NCBI dbSNP database confirmed the build to be GRCh37

# ------------------------------------------------------------------------------
# 1. READ AND PROCESS FILES
# ------------------------------------------------------------------------------
files <- list.files(input_folder, pattern = "Summary_chr", full.names = TRUE)
files <- files[order(as.numeric(gsub("\\D", "", basename(files))))]

results <- list()

message(paste("Found", length(files), "files. Processing..."))

#Checking OR
for(f in files) {
  fname <- basename(f)
  dt <- fread(f, select = c("overall_OR", "serous_hg_OR"))
  
  total_rows <- nrow(dt)
  
  #Count bad rows (OR <= 0)
  bad_overall <- sum(dt$overall_OR <= 0, na.rm = TRUE)
  bad_hgs <- sum(dt$serous_hg_OR <= 0, na.rm = TRUE)
  
  pct_overall <- round((bad_overall / total_rows) * 100, 3)
  pct_hgs <- round((bad_hgs / total_rows) * 100, 3)
  
  results[[fname]] <- data.table(
    File = fname,
    Total_SNPs = total_rows,
    Bad_Overall = bad_overall,
    Pct_Overall = pct_overall,
    Bad_HGS = bad_hgs,
    Pct_HGS = pct_hgs
  )
}

final_report <- rbindlist(results)

print(final_report)

#With ~50% of ORs being <=0 across all files, authors likely meant Beta instead

list_overall <- list()
list_hgs <- list()

for(f in files) {
  message(paste("   Reading:", basename(f)))
  dt <- fread(f)
  
  # ---------------------------------------------------------
  # A. EXTRACT "OVERALL" (All Invasive)
  # ---------------------------------------------------------
  df_ov <- dt[, .(
    snp = `1000G_SNPname`,
    CHR = Chromosome,
    POS = Position,
    ALT = Effect, 
    REF = Baseline,
    BETA = overall_OR,
    SE = overall_SE,
    P = overall_pvalue,
    EAF = EAF
  )]
  list_overall[[f]] <- df_ov
  
  # ---------------------------------------------------------
  # B. EXTRACT "HIGH GRADE SEROUS" (serous_hg)
  # ---------------------------------------------------------
  df_hg <- dt[, .(
    snp = `1000G_SNPname`,
    CHR = Chromosome,
    POS = Position,
    ALT = Effect,
    REF = Baseline,
    BETA = serous_hg_OR, 
    SE = serous_hg_SE,
    P = serous_hg_pvalue,
    EAF = EAF
  )]
  list_hgs[[f]] <- df_hg
}

# ------------------------------------------------------------------------------
# 2. SAVE OUTPUTS
# ------------------------------------------------------------------------------

# --- SAVE OVERALL ---
message("Merging and Saving Overall Dataset...")
final_overall <- rbindlist(list_overall)
fwrite(final_overall, out_overall, sep = "\t")

# --- SAVE HIGH GRADE SEROUS ---
message("Merging and Saving HGS Dataset...")
final_hgs <- rbindlist(list_hgs)
fwrite(final_hgs, out_hgs, sep = "\t")

message("Build 37 saved.")

# ------------------------------------------------------------------------------
# LIFTOVER FUNCTION
# ------------------------------------------------------------------------------

chain_path <- "Data_raw/hg19ToHg38.over.chain"
chain <- import.chain(chain_path)
out_overall_38 <- "Data_raw/Ovarian_Cancer_Overall_Build38.tsv.gz"
out_hgs_38 <- "Data_raw/Ovarian_Cancer_HGS_Build38.tsv.gz"

perform_liftover <- function(dt, chain) {
  chromosomes <- dt$CHR
  if(!any(grepl("chr", head(chromosomes)))) {
    chromosomes <- paste0("chr", chromosomes)
  }
  chromosomes <- gsub("chr23", "chrX", chromosomes)
  chromosomes <- gsub("chr24", "chrY", chromosomes)
  
  gr_37 <- GRanges(
    seqnames = chromosomes,
    ranges = IRanges(start = as.numeric(dt$POS), end = as.numeric(dt$POS)),
    strand = "*" 
  )
  
  message("  Lifting over to hg38...")
  gr_38_list <- liftOver(gr_37, chain)
  n_mapped <- elementNROWS(gr_38_list)
  
  idx_unmapped <- which(n_mapped == 0)
  idx_split <- which(n_mapped > 1) 
  keep_idx <- which(n_mapped == 1)
  
  message("  --- LIFTOVER SUMMARY ---")
  message(paste("  Total Lost:  ", length(idx_unmapped) + length(idx_split)))
  message(paste("  Unmapped (0):", length(idx_unmapped)))
  message(paste("  Split (>1):  ", length(idx_split)))
  message(paste("  Mapped (1):  ", length(keep_idx)))
  
  gr_38 <- unlist(gr_38_list[keep_idx])
  
  dt_clean <- dt[keep_idx, ]
  dt_clean[, POS := start(gr_38)]
  dt_clean[, CHR := gsub("chr", "", as.character(seqnames(gr_38)))]
  
  return(dt_clean)
}

message("Applying Liftover to Overall Dataset...")
final_overall_38 <- perform_liftover(final_overall, chain)

message(paste("Saving Final Build 38:", out_overall_38))
fwrite(final_overall_38, out_overall_38, sep = "\t")

message("Applying Liftover to HGS Dataset...")
final_hgs_38 <- perform_liftover(final_hgs, chain)

message(paste("Saving Final Build 38:", out_hgs_38))
fwrite(final_hgs_38, out_hgs_38, sep = "\t")

message("Done.")
