library(data.table)

target_snps <- c("rs7605725", "rs6062498")
raw_file <- "Data_raw/28067908-GCST004133-EFO_0000729.h.tsv.gz"
cols <- names(fread(raw_file, nrows=0))
print(cols)

raw_cols_to_keep <- c("hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "effect_allele", "other_allele", "beta")

dt_raw <- fread(raw_file, select = raw_cols_to_keep)
dt_filtered <- dt_raw[hm_rsid %in% target_snps]

print(dt_filtered)