library(data.table)

target_snp <- c("rs7017073", "rs7650602")


raw_file  <- "Data_Raw/oncoarray_bcac_public_release_oct17.txt.gz"
conv_file <- "Data_Raw/BC_ERpos_Combined_Build37.tsv.gz"

raw_cols_to_keep <- c("phase3_1kg_id", "a0", "a1", "bcac_onco_icogs_gwas_erpos_beta")
dt_raw <- fread(raw_file, select = raw_cols_to_keep)[phase3_1kg_id %like% target_snp]

setnames(dt_raw, 
         old = c("phase3_1kg_id", "a0", "a1", "bcac_onco_icogs_gwas_erpos_beta"), 
         new = c("SNP_RAW", "RAW_REF(a0)", "RAW_ALT(a1)", "RAW_BETA"))

message("-> Found in RAW:")
print(dt_raw)

conv_cols_to_keep <- c("snp", "REF", "ALT", "BETA")
dt_conv <- fread(conv_file, select = conv_cols_to_keep)[snp %like% target_snp]

message("-> Found in CONVERTED:")
print(dt_conv)

setnames(dt_conv, 
         old = c("snp", "REF", "ALT", "BETA"), 
         new = c("SNP", "CONV_REF", "CONV_ALT", "CONV_BETA"))

comparison_df <- cbind(dt_raw, dt_conv)
print(comparison_df)
