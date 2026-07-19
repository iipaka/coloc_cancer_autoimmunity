library(data.table)

target_snps <- "rs7513707"
raw_file <- "Data_raw/24390342-GCST002318-EFO_0000685.h.tsv.gz"
cols <- names(fread(raw_file, nrows=0))
print(cols)

raw_cols_to_keep <- c("hm_rsid", "hm_effect_allele", "hm_other_allele", "beta", "hm_odds_ratio")

dt_raw <- fread(raw_file, select = raw_cols_to_keep)
dt_filtered <- dt_raw[hm_rsid %in% target_snps]

print(dt_filtered)