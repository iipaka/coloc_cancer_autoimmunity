# ==============================================================================
# Purpose: Generate diagnostic plots (P-value & Z-score); GWAS-GWAS
# ==============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

out_dir <- "Output/Diagnostic_Plots_GWAS_precoloc"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gwas_info <- read.csv("Data_raw/gwas_summarystats.csv")
p_threshold <- "1e-06"

lead_coloc_snps <- read.csv(paste0("Output/sig_snps_threshold_", p_threshold, ".csv"))

plot_theme <- theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 18),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )

# ==============================================================================
# LOOP THROUGH THE HITS
# ==============================================================================

for(i in 1:nrow(lead_coloc_snps)) {
  
  rsid <- lead_coloc_snps$rsID[i]
  trait1_name <- lead_coloc_snps$Cancer_Type[i] 
  trait2_name <- lead_coloc_snps$Name[i]

  message(sprintf("[%d/%d] Plotting diagnostics for %s (%s vs %s)", i, nrow(lead_coloc_snps), rsid, trait1_name, trait2_name))  
  
  trait1_gcst <- gwas_info$name[gwas_info$Disease == trait1_name]
  trait2_gcst <- gwas_info$name[gwas_info$Disease == trait2_name]
  
  file_trait1 <- file.path("Data_processed", paste0(rsid, "_", trait1_gcst, "_", p_threshold, ".rds"))
  file_trait2 <- file.path("Data_processed", paste0(rsid, "_", trait2_gcst, "_", p_threshold, ".rds"))
  
  d_trait1 <- readRDS(file_trait1)
  d_trait2 <- readRDS(file_trait2)
  
  df_trait1 <- data.frame(
    snp = d_trait1$snp,
    rsid = d_trait1$rsid,
    p_t1 = d_trait1$pvalues,
    z_t1 = d_trait1$beta / sqrt(d_trait1$varbeta)
  )
  
  df_trait2 <- data.frame(
    snp = d_trait2$snp,
    rsid = d_trait2$rsid,
    p_t2 = d_trait2$pvalues,
    z_t2 = d_trait2$beta / sqrt(d_trait2$varbeta)
  )
  
  df_plot <- inner_join(df_trait1, df_trait2, by = "snp")
  
  #Calculate the -log10 P-value (adding a small constant to prevent log(0) errors)
  df_plot$logp_t1 <- -log10(df_plot$p_t1 + 1e-300)
  df_plot$logp_t2 <- -log10(df_plot$p_t2 + 1e-300)
  
  lead_snp_data <- df_plot[grepl(rsid, df_plot$rsid.x), ]
  
  p1 <- ggplot(df_plot, aes(x = logp_t1, y = logp_t2)) +
    geom_point(alpha = 0.6, color = "darkblue", size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", alpha = 0.5) +
    labs(title = "P-value Correlation",
         x = paste0("-log10(P) [", trait1_name, "]"),
         y = paste0("-log10(P) [", trait2_name, "]")) +
    plot_theme
  

  p2 <- ggplot(df_plot, aes(x = z_t1, y = z_t2)) +
    geom_point(alpha = 0.6, color = "darkgreen", size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", alpha = 0.5) +
    labs(title = "Z-score (Effect Direction)",
         x = paste0("Z-score [", trait1_name, "]"),
         y = paste0("Z-score [", trait2_name, "]")) +
    plot_theme
  
  combined_plot <- p1 + p2 + 
    plot_annotation(title = paste0("Diagnostic Plots: ", rsid, " (", trait1_name, " vs ", trait2_name, ")"),
                    theme = theme(plot.title = element_text(size = 20)))
  
  safe_rsid <- gsub(":", "_", rsid)
  plot_file <- file.path(out_dir, paste0("Diagnostic_", safe_rsid, "_", trait1_name, "_vs_", trait2_name, ".png"))
  
  ggsave(plot_file, combined_plot, width = 10, height = 5, dpi = 300)
  message("  -> Saved: ", plot_file)
}

message("\nAll diagnostic plots complete.")