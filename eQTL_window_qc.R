library(data.table)

chr <- "12"
pos <- 111446804 
window_start <- pos - 1000000/2 
window_end <- pos + 1000000/2

# Query the raw eQTLGen file using tabix
# (Make sure to replace the filename with your actual raw eQTLGen file)
raw_eqtl_file <- "/rds/user/xz418/hpc-work/Data_eQTLGen/Blood_eQTLGen_hg38_sorted.tsv.gz"
region_query <- paste0(chr, ":", window_start, "-", window_end)

cmd <- paste("tabix", raw_eqtl_file, region_query)

qc_data <- fread(cmd)

colnames(qc_data) <- c("CHR38", "BP38", "SNPID", "CHR19", "BP19", "REF", "ALT", "ALT_FREQ", 
                       "P", "Zscore", "Gene_id", "Gene_name", "GeneChr", "GenePos", 
                       "NrCohorts", "NrSamples", "FDR", "BonferroniP", "BETA", "SE")

unique_genes <- unique(qc_data$Gene_name)

cat("Found", length(unique_genes), "unique genes tested in this window:\n")
print(unique_genes)