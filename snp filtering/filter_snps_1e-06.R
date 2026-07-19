# ==========================================
# Purpose: Filters for SNPs significant in cancer and immune-mediated diseases
# Input: Excel file in "Data_raw"
# Output: CSV file in "Output"
# ==========================================

library(readxl)
library(dplyr)
library(purrr)
library(rtracklayer)
library(GenomicRanges)

file_path <- "Data_raw/genes-3623764-supplementary.xlsx"
chain_path <- "Data_raw/hg19ToHg38.over.chain"
target_sheets <- c("Table S3", "TableS4", "Table S5", "Table S6", "Table S7", "Table S8", "Table S9")

p_threshold<-1e-6

#Function to process excel sheet
process_gwas_sheet<-function(sheet_name){
  df<-read_excel(file_path, sheet=sheet_name, skip=2)
  colnames(df)[1]<-"rsID"
  colnames(df)[2]<-"chr"
  colnames(df)[3]<-"pos"
  colnames(df)[9]<-"Cancer_Pval"
  colnames(df)[13]<-"Immune_Pval"
  filtered_df<-df%>%
    filter(Cancer_Pval<p_threshold & Immune_Pval<p_threshold) %>%
    mutate(
      Source_Sheet = sheet_name,
      Cancer_Type = case_when(
        sheet_name == "Table S3" ~ "BC",
        sheet_name == "TableS4" ~ "BC_erpos",
        sheet_name == "Table S5" ~ "BC_erneg",
        sheet_name == "Table S6" ~ "OC",
        sheet_name == "Table S7" ~ "OC_HGS",
        sheet_name == "Table S8" ~ "PC",
        sheet_name == "Table S9" ~ "EC")
    )
  return(filtered_df)
}


all_significant_snps<-map_dfr(target_sheets, process_gwas_sheet)

chain <- import.chain(chain_path)

chromosomes <- paste0("chr", all_significant_snps$chr)

gr_37 <- GRanges(
  seqnames = chromosomes,
  ranges = IRanges(start = as.numeric(all_significant_snps$pos), 
                   end = as.numeric(all_significant_snps$pos)),
  strand = "*"
)

gr_38_list <- liftOver(gr_37, chain)
n_mapped <- elementNROWS(gr_38_list)

idx_unmapped <- which(n_mapped == 0)
idx_split <- which(n_mapped > 1) 
keep_idx <- which(n_mapped == 1)

message("--- LIFTOVER SUMMARY ---")
message(paste("Total Lost:  ", length(idx_unmapped) + length(idx_split)))
message(paste("Unmapped (0):", length(idx_unmapped)))
message(paste("Split (>1):  ", length(idx_split)))
message(paste("Mapped (1):  ", length(keep_idx)))

gr_38 <- unlist(gr_38_list[keep_idx])

final_snps <- all_significant_snps[keep_idx, ]
final_snps$pos <- start(gr_38)
final_snps$chr <- gsub("chr", "", as.character(seqnames(gr_38)))

output_filename <- paste0("Output/sig_snps_threshold_", p_threshold, ".csv")
write.csv(final_snps, output_filename, row.names=FALSE)
message("Extraction complete.")
