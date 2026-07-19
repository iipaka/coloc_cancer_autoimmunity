# ==============================================================================
# Purpose: Assign nearest protein-coding gene to colocalisation hits
# ==============================================================================

library(dplyr)
library(biomaRt)

coloc_res <- read.csv("Output/coloc_results_summary_1e-06.csv")

lead_coloc_snps <- coloc_res %>% filter(PP_H4 > 0.75)

message("Found ", nrow(lead_coloc_snps), " lead SNPs after coloc.")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror="useast")

get_nearest_gene <- function(chr, pos) {

  window <- 1000000
  start_pos <- pos - window
  end_pos <- pos + window
  
  genes <- getBM(
    attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position', 'gene_biotype'),
    filters = c('chromosome_name', 'start', 'end', 'biotype'),
    values = list(chr, start_pos, end_pos, 'protein_coding'),
    mart = ensembl
  )
  
  genes <- genes %>% filter(hgnc_symbol != "")
  
  if (nrow(genes) == 0) return(NA)
  
  genes <- genes %>%
    rowwise() %>%
    mutate(
      dist_to_gene = case_when(
        pos >= start_position & pos <= end_position ~ 0,
        pos < start_position ~ start_position - pos,
        pos > end_position ~ pos - end_position
      )
    ) %>%
    ungroup() %>%
    arrange(dist_to_gene)
  
  return(genes$hgnc_symbol[1])
}

lead_coloc_snps$Nearest_gene <- mapply(get_nearest_gene, lead_coloc_snps$CHR, lead_coloc_snps$POS)

final_table <- lead_coloc_snps %>% 
  dplyr::select(SNP, CHR, POS, Cancer_trait, Immune_trait, PP_H4, Nearest_gene)

print(final_table)

write.csv(final_table, "Output/coloc_lead_snps_nearest_gene_1e-06.csv", row.names = FALSE)
