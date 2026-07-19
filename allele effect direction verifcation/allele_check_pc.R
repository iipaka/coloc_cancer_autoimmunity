library(data.table)

target_snps <- c("rs7605725", "rs6062498", "rs2296763")
raw_file <- "Data_raw/29892016-GCST006085-EFO_0001663.h.tsv.gz"
raw_cols_to_keep <- c("hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "effect_allele", "other_allele", "beta")

dt_raw <- fread(raw_file, select = raw_cols_to_keep)
dt_filtered <- dt_raw[snp %like% target_snps]

print(dt_filtered)
