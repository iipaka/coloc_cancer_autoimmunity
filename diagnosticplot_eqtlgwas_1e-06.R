# ==============================================================================
# Purpose: Generate LocusCompare diagnostic plots (P-value & Z-score scatters)
# ==============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

# 1. Setup
out_dir <- "Output/Diagnostic_Plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gwas_info <- read.csv("Data_raw/gwas_summarystats.csv")
p_threshold <- "1e-06"
coloc_res <- read.csv(paste0("Output/coloc_results_eqtl_gwas_summary_", p_threshold, ".csv"))
lead_coloc_snps <- coloc_res %>% filter(PP_H4 > 0.75)

# 2. Plotting Themes
plot_theme <- theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11, face = "bold"))

# ==============================================================================
# LOOP THROUGH THE HITS
# ==============================================================================

for(i in 1:nrow(lead_coloc_snps)) {
  
  rsid <- lead_coloc_snps$Target_SNP[i]
  trait_name <- lead_coloc_snps$GWAS_Trait[i]
  gene_name <- lead_coloc_snps$eQTL_Gene[i]
  pp4 <- round(lead_coloc_snps$PP_H4[i], 3)
  
  message(sprintf("[%d/%d] Plotting diagnostics for %s (%s vs %s)", i, nrow(lead_coloc_snps), rsid, trait_name, gene_name))
  
  # A. Find correct files
  trait_gcst <- gwas_info$name[gwas_info$Disease == trait_name]
  file_gwas <- file.path("Data_processed", paste0(rsid, "_", trait_gcst, "_", p_threshold, ".rds"))
  file_eqtl <- file.path("Data_processed", paste0("eqtl_", rsid, "_", gene_name, ".rds"))
  
  
  # B. Load data
  d_gwas <- readRDS(file_gwas)
  d_eqtl <- readRDS(file_eqtl)
  
  # C. Extract and calculate Z-scores (Z = beta / standard_error)
  # coloc lists store varbeta (SE^2), so SE is sqrt(varbeta)
  df_gwas <- data.frame(
    snp = d_gwas$snp,
    p_gwas = d_gwas$pvalues,
    z_gwas = d_gwas$beta / sqrt(d_gwas$varbeta)
  )
  
  df_eqtl <- data.frame(
    snp = d_eqtl$snp,
    p_eqtl = d_eqtl$pvalues,
    z_eqtl = d_eqtl$beta / sqrt(d_eqtl$varbeta)
  )
  
  # D. Merge on exact SNP matches
  df_plot <- inner_join(df_gwas, df_eqtl, by = "snp")

  
  # Calculate -log10 P-values
  df_plot$logp_gwas <- -log10(df_plot$p_gwas + 1e-300)
  df_plot$logp_eqtl <- -log10(df_plot$p_eqtl + 1e-300)
  
  # E. Generate P-Value Plot
  p1 <- ggplot(df_plot, aes(x = logp_gwas, y = logp_eqtl)) +
    geom_point(alpha = 0.6, color = "darkblue", size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", alpha = 0.5) +
    labs(title = "P-value Correlation",
         subtitle = paste0("Coloc PP.H4: ", pp4),
         x = paste0("-log10(P) [", trait_name, "]"),
         y = paste0("-log10(P) [", gene_name, " eQTL]")) +
    plot_theme
  
  # F. Generate Z-Score Plot
  p2 <- ggplot(df_plot, aes(x = z_gwas, y = z_eqtl)) +
    geom_point(alpha = 0.6, color = "darkgreen", size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", alpha = 0.5) +
    labs(title = "Z-score (Effect Direction)",
         x = paste0("Z-score [", trait_name, "]"),
         y = paste0("Z-score [", gene_name, " eQTL]")) +
    plot_theme
  
  # G. Combine and Save
  combined_plot <- p1 + p2 + 
    plot_annotation(title = paste0("Diagnostic Plots: ", rsid, " (", trait_name, " vs ", gene_name, ")"))
  
  safe_rsid <- gsub(":", "_", rsid)
  plot_file <- file.path(out_dir, paste0("Diagnostic_", safe_rsid, "_", trait_name, "_vs_", gene_name, ".png"))
  
  ggsave(plot_file, combined_plot, width = 10, height = 5, dpi = 300)
  message("  -> Saved: ", plot_file)
}

message("\nAll diagnostic plots complete.")