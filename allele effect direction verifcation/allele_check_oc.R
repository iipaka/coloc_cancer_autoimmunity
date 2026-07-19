library(data.table)

target_snps <- "rs10089868"
raw_file <- "Data_raw/Ovarian_Cancer_Overall_Build38.tsv.gz"
cols <- names(fread(raw_file, nrows=0))
print(cols)

raw_cols_to_keep <- c("snp", "ALT", "REF", "BETA", "POS")

dt_raw <- fread(raw_file, select = raw_cols_to_keep)
dt_filtered <- dt_raw[snp %like% target_snps]

print(dt_filtered)
