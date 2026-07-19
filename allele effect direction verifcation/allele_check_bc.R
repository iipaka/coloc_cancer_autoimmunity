library(data.table)

target_snps <- "rs7017073"
raw_file <- "Data_raw/29059683-GCST004988-EFO_0000305.h.tsv.gz"
cols <- names(fread(raw_file, nrows=0))
print(cols)

raw_cols_to_keep <- c("hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta")

dt_raw <- fread(raw_file, select = raw_cols_to_keep)
dt_filtered <- dt_raw[hm_rsid %in% target_snps]

print(dt_filtered)
