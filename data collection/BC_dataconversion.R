# ==============================================================================
# Purpose: Standardise Breast Cancer GWAS
# ==============================================================================

library(data.table)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(R.utils)

# ------------------------------------------------------------------------------
# 1. SETUP
# ------------------------------------------------------------------------------

input_file <- "Data_Raw/oncoarray_bcac_public_release_oct17.txt.gz"
dt_preview <- fread(input_file, nrows = 5) #Read first 5 rows to inspect structure

message("\n--- COLUMN NAMES ---")
print(names(dt_preview))

message("\n--- DATA PREVIEW (First 5 Rows) ---")
print(dt_preview)

#Chr: col 3, Pos: col 4

out_erpos <- "Data_Raw/BC_ERpos_Combined_Build37.tsv.gz"
out_erneg <- "Data_Raw/BC_ERneg_Combined_Build37.tsv.gz"

message(paste("Reading master file:", input_file))
dt <- fread(input_file)

# ------------------------------------------------------------------------------
# 2. DEFINE COLUMN MAPPING
# ------------------------------------------------------------------------------

# --- ER POSITIVE  ---
message("Processing ER+...")

#Create new table for ER+ data by selecting specific columns from the master table and renaming them for standardisation
df_erpos <- dt[, .(
  snp = phase3_1kg_id,
  CHR = chr,
  POS = position_b37,
  REF = a0,
  ALT = a1,
  EAF = bcac_onco_icogs_gwas_erpos_eaf_controls,
  BETA = bcac_onco_icogs_gwas_erpos_beta,
  SE = bcac_onco_icogs_gwas_erpos_se,
  P = bcac_onco_icogs_gwas_erpos_P1df
)]

n_before <- nrow(df_erpos)

df_erpos <- df_erpos[!is.na(BETA) & !is.na(SE)] #Filter variants missing BETA or SE
fwrite(df_erpos, out_erpos, sep = "\t")

n_after <- nrow(df_erpos)
n_lost <- n_before - n_after
percent_lost <- (n_lost / n_before) * 100

message(paste0("QC removed ", n_lost, "variants (", round(percent_lost, 2), "% of total)."))

# --- ER NEGATIVE ---
message("Processing ER-...")
df_erneg <- dt[, .(
  snp = phase3_1kg_id,
  CHR = chr,
  POS = position_b37,
  REF = a0,
  ALT = a1,
  EAF = bcac_onco_icogs_gwas_erneg_eaf_controls,
  BETA = bcac_onco_icogs_gwas_erneg_beta,
  SE = bcac_onco_icogs_gwas_erneg_se,
  P = bcac_onco_icogs_gwas_erneg_P1df
)]

n_before <- nrow(df_erneg)

df_erneg <- df_erneg[!is.na(BETA) & !is.na(SE)]
fwrite(df_erneg, out_erneg, sep = "\t")

n_after <- nrow(df_erneg)
n_lost <- n_before - n_after
percent_lost <- (n_lost / n_before) * 100

message(paste0("QC removed ", n_lost, "variants (", round(percent_lost, 2), "% of total)."))

# ------------------------------------------------------------------------------
# 3. GENOME BUILD MAPPING
# ------------------------------------------------------------------------------

files_to_process <- c(
  "Data_Raw/BC_ERpos_Combined_Build37.tsv.gz",
  "Data_Raw/BC_ERneg_Combined_Build37.tsv.gz")

suffix <- "_Build38.tsv.gz"

chain_url <- "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
chain_file <- "Data_Raw/hg19ToHg38.over.chain.gz"

if (!file.exists(chain_file)) {
  message("Downloading Chain File (hg19 -> hg38)...")
  download.file(chain_url, chain_file)
}

message("Importing Chain file...")
#gunzip("Data_Raw/hg19ToHg38.over.chain.gz", remove=FALSE) #only needed the first time
chain <- import.chain("Data_Raw/hg19ToHg38.over.chain")

# ------------------------------------------------------------------------------
# 4. LIFTOVER
# ------------------------------------------------------------------------------

#For each of the build 37 files:
for (f_path in files_to_process) {
  
  message(paste0("\nProcessing: ", basename(f_path)))
  
  dt <- fread(f_path) #Read file
  
  chromosomes <- dt$CHR #Extract chr column
  
  #Check if chr names start with "chr", if not, add "chr" (needed by liftover)
  if(!any(grepl("chr", head(chromosomes)))) {
    chromosomes <- paste0("chr", chromosomes)
  }
  
  #Convert chr23 and chr24 to chrX and chrY respectively
  chromosomes <- gsub("chr23", "chrX", chromosomes)
  chromosomes <- gsub("chr24", "chrY", chromosomes)
  
  #Construct GenomicRanges object representing the original build 37 coordinates for each variant
  gr_37 <- GRanges(
    seqnames = chromosomes,
    ranges = IRanges(start = dt$POS, end = dt$POS),
    strand = "*" 
  )
  
  message("   Lifting over to hg38...")
  gr_38_list <- liftOver(gr_37, chain)
  
  n_mapped <- elementNROWS(gr_38_list) #Count how many hg38 locations each original variant mapped to (0=unmapped, 1=successful, >1=split mapping).
  
  idx_unmapped <- which(n_mapped == 0)
  idx_split <- which(n_mapped > 1) 
  keep_idx <- which(n_mapped == 1)
  
  message("\n--- LIFTOVER SUMMARY ---")
  message(paste("Total Lost:      ", length(idx_unmapped) + length(idx_split)))
  message(paste("Unmapped (0): ", length(idx_unmapped), "(Deleted/Gapped in hg38)"))
  message(paste("Split (>1):   ", length(idx_split),    "(Maps to multiple locations)"))
  message(paste0("   Mapped ", length(keep_idx), " variants. Lost ", nrow(dt) - length(keep_idx), "."))
  
  gr_38 <- unlist(gr_38_list[keep_idx]) #Extract genomic coordinates for variants that mapped successfully

  dt_clean <- dt[keep_idx, ] #Subset original data table to keep only successfully mapped variants
  dt_clean[, POS := start(gr_38)] #Replace old POS column with new hg38 pos
  dt_clean[, CHR := gsub("chr", "", as.character(seqnames(gr_38)))] #Strip chr prefix
  
  new_filename <- gsub("_Build38.tsv.gz", "", f_path) 
  new_filename <- gsub("_Build37.tsv.gz", "_Build38.tsv.gz", new_filename)
  
  message(paste("   Saving to:", new_filename))
  fwrite(dt_clean, new_filename, sep = "\t")
}

message("Done.")
