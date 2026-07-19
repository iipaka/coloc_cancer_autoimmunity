# ==============================================================================
# Purpose: Generate narrow-window diagnostic plots (P-value & Z-score); GWAS-EQTL
# ==============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

p_threshold <- "1e-06"
out_dir <- paste0("Output/Diagnostic_Plots_eQTL_", p_threshold, "_NarrowWindow")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_plots <- data.frame(
  rsid = "rs3184504",
  trait_gcst = "GCST006464",
  trait_name = "EC", 
  gene_name = "PPTC7",
  pos_min = 111300000, 
  pos_max = 111600000, 
  stringsAsFactors = FALSE
)

plot_theme <- theme_bw(base_family = "Arial") + 
  theme(
    text = element_text(family = "Arial"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 11, face = "bold", family = "Arial")
  )

# ==============================================================================
# RUN THE TARGETED PLOT
# ==============================================================================

for(i in 1:nrow(target_plots)) {
  
  rsid <- target_plots$rsid[i]
  trait_gcst <- target_plots$trait_gcst[i]
  trait_name <- target_plots$trait_name[i]
  gene_name <- target_plots$gene_name[i]
  pos_min <- target_plots$pos_min[i]
  pos_max <- target_plots$pos_max[i]
  
  message(sprintf("Plotting narrow window for %s (%s vs %s)", rsid, trait_name, gene_name))
  
  file_gwas <- file.path("Data_processed", paste0(rsid, "_", trait_gcst, "_", p_threshold, "_centered.rds"))
  file_eqtl <- file.path("Data_processed", paste0("eqtl_", rsid, "_", gene_name, "_centered.rds"))
  
  d_gwas <- readRDS(file_gwas)
  d_eqtl <- readRDS(file_eqtl)
  
  df_gwas <- data.frame(
    snp = d_gwas$snp,
    position = d_gwas$position,
    p_gwas = d_gwas$pvalues,
    z_gwas = d_gwas$beta / sqrt(d_gwas$varbeta)
  )
  
  df_eqtl <- data.frame(
    snp = d_eqtl$snp,
    p_eqtl = d_eqtl$pvalues,
    z_eqtl = d_eqtl$beta / sqrt(d_eqtl$varbeta)
  )
  
  df_plot <- inner_join(df_gwas, df_eqtl, by = "snp") %>%
    filter(position >= pos_min & position <= pos_max)
  
  df_plot$logp_gwas <- -log10(df_plot$p_gwas + 1e-300)
  df_plot$logp_eqtl <- -log10(df_plot$p_eqtl + 1e-300)
  
  p1 <- ggplot(df_plot, aes(x = logp_gwas, y = logp_eqtl)) +
    geom_point(alpha = 0.6, color = "darkblue", size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", alpha = 0.5) +
    labs(x = paste0("-log10(P) [", trait_name, "]"),
         y = paste0("-log10(P) [", gene_name, " eQTL]")) +
    plot_theme
  
  p2 <- ggplot(df_plot, aes(x = z_gwas, y = z_eqtl)) +
    geom_point(alpha = 0.6, color = "darkgreen", size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", alpha = 0.5) +
    labs(x = paste0("Z-score [", trait_name, "]"),
         y = paste0("Z-score [", gene_name, " eQTL]")) +
    plot_theme
  
  combined_plot <- p1 + p2 & 
    theme(text = element_text(family = "Arial"))
  
  safe_rsid <- gsub(":", "-", rsid)
  safe_gene <- gsub(":", "-", gene_name)
  plot_file <- file.path(out_dir, paste0("NarrowDiag_", safe_rsid, "_", trait_name, "_vs_", safe_gene, ".png"))
  
  ggsave(plot_file, combined_plot, width = 10, height = 5, dpi = 300)
  message("  -> Saved: ", plot_file)
}

message("\nAll narrow diagnostic plots complete.")
