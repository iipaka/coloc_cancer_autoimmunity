# ==============================================================================
# Purpose: Generate Manhattan plots (GWAS <> eQTL) on Centered Data (POST-COLOC)
# ==============================================================================
library(ggplot2)
library(dplyr)
library(data.table)
library(patchwork)

p_threshold <- "1e-06"

out_dir <- paste0("Output/Plots_eQTL_", p_threshold, "_centered_postcoloc")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sig_snps <- read.csv("Output/recentered_targets.csv")
gwas_info <- read.csv("Data_raw/gwas_summarystats.csv")

coloc_res <- read.csv(paste0("Output/coloc_results_eqtl_gwas_summary_", p_threshold, "_centered.csv"))
lead_coloc_snps <- coloc_res %>% filter(PP_H4 > 0.75)

# ==============================================================================
# PLOTTING FUNCTION
# ==============================================================================

manhattan_theme <- theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 10, face = "plain")
  )

gwas_threshold <- -log10(as.numeric(p_threshold))

manhattan_plot_fx <- function(data, title, shared_pos, target_pos, dot_colour, x_min, x_max) {
  
  ggplot(data, aes(x = POS / 1e6, y = logP)) +

    geom_vline(xintercept = target_pos / 1e6, linetype = "dashed", alpha = 0.5) +
    
    geom_point(data = filter(data, POS %in% shared_pos), 
               aes(color = "Common to both"), alpha = 0.7, size = 1.75, stroke = 0.5) +
    
    geom_point(data = filter(data, !(POS %in% shared_pos)), 
               aes(color = Dataset), alpha = 0.6, size = 1.75, stroke = 0.5) +
    
    geom_hline(yintercept = gwas_threshold, linetype = "dotted", color = "black") +
    annotate("text", x = Inf, y = gwas_threshold, label = "p-value threshold", 
             vjust = -0.5, hjust = 1.05, color = "black", size = 2.5) +
    
    scale_color_manual(
      name = "SNP", 
      breaks = c(title, "Common to both"), 
      values = setNames(c("pink", dot_colour), c("Common to both", title))
    ) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    labs(title = title, x = expression("BP38 / 10"^6), y = expression("-log"[10]*"(P-value)")) +
    manhattan_theme
}

# ==============================================================================
# LOOP THROUGH SNPs
# ==============================================================================

for(i in 1:nrow(sig_snps)) {
  
  target_snp <- sig_snps$New_Centre_SNP[i]
  clean_snp <- sub(":.*", "", target_snp)
  target_chr <- sig_snps$CHR[i]
  target_pos <- sig_snps$New_Centre_POS[i]
  trait_name <- sig_snps$Stronger_Trait[i]
  
  region_hits <- lead_coloc_snps %>% 
    filter(Target_SNP == target_snp & GWAS_Trait == trait_name)
  
  if(nrow(region_hits) == 0) {
    message(paste0("\n[", i, "/", nrow(sig_snps), "] Skipping ", clean_snp, " (No hits > 0.75)"))
    next
  }
  
  trait_gcst <- gwas_info$name[gwas_info$Disease == trait_name]
  
  file_gwas <- paste0("Data_processed/", target_snp, "_", trait_gcst, "_", p_threshold, "_centered.rds")
  data_gwas <- readRDS(file_gwas)
  
  df_gwas <- data.frame(
    POS = data_gwas$position,
    P = data_gwas$pvalues,
    Dataset = trait_name 
  )
  
  for(j in 1:nrow(region_hits)) {
    
    gene_name <- region_hits$eQTL_Gene[j]
    pp4       <- round(region_hits$PP_H4[j], 3)
    
    message(paste0("\n  -> Plotting hit: ", trait_name, " vs ", gene_name, " eQTL (PP.H4: ", pp4, ")"))
    
    file_eqtl <- paste0("Data_processed/eqtl_", target_snp, "_", gene_name, "_centered.rds")
    data_eqtl <- readRDS(file_eqtl)
    
    df_eqtl <- data.frame(
      POS = data_eqtl$position,
      P = data_eqtl$pvalues,
      Dataset = paste0(gene_name, " eQTL")
    )
    
    shared_positions <- intersect(df_gwas$POS, df_eqtl$POS)
    
    df_plot <- rbind(df_gwas, df_eqtl)
    df_plot$logP <- -log10(df_plot$P + 1e-300) 
    
    window_min <- min(df_plot$POS, na.rm = TRUE) / 1e6
    window_max <- max(df_plot$POS, na.rm = TRUE) / 1e6
    
    p1 <- manhattan_plot_fx(filter(df_plot, Dataset == trait_name), 
                            trait_name, shared_positions, target_pos, "red", window_min, window_max)
    
    p2 <- manhattan_plot_fx(filter(df_plot, Dataset == paste0(gene_name, " eQTL")), 
                            paste0(gene_name, " eQTL"), shared_positions, target_pos, "blue", window_min, window_max)
    
    combined_plot <- (p1 / p2) + 
      plot_annotation(
        title = paste0("Manhattan Plot: ", clean_snp, " (Chromosome ", target_chr, ")"),
        subtitle = paste0("PP_H4: ", pp4, " | eQTL Gene: ", gene_name)
      )
    
    plot_filename <- file.path(out_dir, paste0("Plot_", clean_snp, "_", trait_name, "_vs_", gene_name, "_eQTL.png"))
    
    ggsave(plot_filename, plot = combined_plot, width = 10, height = 8, dpi = 300)
    message(paste("  -> Saved:", plot_filename))
  }
}