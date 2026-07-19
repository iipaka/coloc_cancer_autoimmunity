# ==============================================================================
# Purpose: Merge lead coloc hits with original p-values to determine the 
#          stronger trait for eQTL analysis.
# ==============================================================================

library(dplyr)

lead_snps <- read.csv("Output/coloc_lead_snps_nearest_gene_1e-06.csv")
sig_info <- read.csv("Output/sig_snps_threshold_1e-06.csv") #Master list of sig snps to retrieve disease p-values

#Extract and rename specific columns from master sig snps list to match lead snps table from nearest gene analysis
sig_pvals <- sig_info %>%
  dplyr::select(
    SNP = rsID,
    Cancer_trait = Cancer_Type,
    Immune_trait = Name,
    Cancer_P = Cancer_Pval,
    Immune_P = Immune_Pval
  )

#Merge original disease p values by joining using snp, cancer and disease name
merged_targets <- lead_snps %>%
  left_join(sig_pvals, by = c("SNP", "Cancer_trait", "Immune_trait"))

base_targets <- merged_targets %>%
  rowwise() %>%
  mutate(
    #Create col to check if cancer p-value is smaller. If true, assign cancer name, else immune name
    Stronger_Trait = ifelse(Cancer_P < Immune_P, Cancer_trait, Immune_trait),
    #Create col to keep smaller p-value
    Stronger_Pval = min(Cancer_P, Immune_P)
  ) %>%
  ungroup()

t1d_additional <- base_targets %>%
  filter(Stronger_Trait == "T1D") %>%
  mutate(
    Stronger_Trait = Cancer_trait,
    Stronger_Pval = Cancer_P # Ensure the p-value updates to the cancer p-value for downstream eQTL steps
  )

final_targets <- bind_rows(base_targets, t1d_additional)

write.csv(final_targets, "Output/eQTL_coloc_gwas_targets_1e-06.csv", row.names = FALSE)