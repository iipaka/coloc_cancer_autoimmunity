library(data.table)

target_snps <- c("rs7650602","rs938648", "rs10089868")

raw_file <- "Data_raw/GCST90446792.h.tsv.gz"
cols <- names(fread(raw_file, nrows=0))
print(cols)

raw_cols_to_keep <- c("rsid", "effect_allele", "other_allele", "beta")

dt_raw <- fread(raw_file, select = raw_cols_to_keep)
dt_filtered <- dt_raw[rsid %in% target_snps]

print(dt_filtered)