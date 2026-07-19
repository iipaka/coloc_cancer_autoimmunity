library(data.table)

target_snps <- c("rs7650602","rs10089868", "rs938648")
raw_file <- "Data_raw/28067908-GCST004132-EFO_0000384.h.tsv.gz"
cols <- names(fread(raw_file, nrows=0))
print(cols)

raw_cols_to_keep <- c("hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "effect_allele", "other_allele", "beta")

dt_raw <- fread(raw_file, select = raw_cols_to_keep)
dt_filtered <- dt_raw[hm_rsid %in% target_snps]

print(dt_filtered)