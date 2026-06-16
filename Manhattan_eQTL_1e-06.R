# ==============================================================================
# Purpose: Dynamically generate regional Manhattan plots (GWAS <> eQTL)
# ==============================================================================
library(ggplot2)
library(dplyr)
library(data.table)
library(patchwork)

p_threshold <- 1e-06

out_dir <- paste0("Output/Plots_eQTL_", p_threshold)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sig_snps <- read.csv("Output/eQTL_coloc_gwas_targets_1e-06.csv")
gwas_info <- read.csv("Data_raw/gwas_summarystats.csv")

coloc_res <- read.csv(paste0("Output/coloc_results_eqtl_gwas_summary_", p_threshold, ".csv"))
lead_coloc_snps <- coloc_res %>% filter(PP_H4 > 0.75)

# ==============================================================================
# PLOTTING
# ==============================================================================

manhattan_theme <- theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 10, face = "plain")
  )

gwas_threshold <- -log10(p_threshold)

manhattan_plot_fx <- function(data, title, shared_pos, target_pos, target_name, dot_colour, x_min, x_max) {
  
  ggplot(data, aes(x = POS / 1e6, y = logP)) +
    
    # Vertical line for target
    geom_vline(xintercept = target_pos / 1e6, linetype = "dashed", alpha = 0.5) +
    
    # Shared points
    geom_point(data = filter(data, POS %in% shared_pos), 
               aes(color = "Common to both"), alpha = 0.7, size = 1.75, stroke = 0.5) +
    
    # Base points
    geom_point(data = filter(data, !(POS %in% shared_pos)), 
               aes(color = Dataset), alpha = 0.6, size = 1.75, stroke = 0.5) +
    
    # Target SNP highlighting
    geom_point(data = filter(data, POS == target_pos), 
               shape = 1, size = 4, color = "black", stroke = 1) +
    
    # Target SNP text
    geom_text(data = filter(data, POS == target_pos), 
              aes(label = target_name), 
              hjust = -0.2, vjust = 0.5, fontface = "bold", size = 3) +
    
    # Significance line and text
    geom_hline(yintercept = gwas_threshold, linetype = "dotted", color = "black") +
    annotate("text", x = Inf, y = gwas_threshold, label = "p-value threshold", 
             vjust = -0.5, hjust = 1.05, color = "black", size = 2.5) +
    
    # Scale and Labels
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
  
  rsid        <- sig_snps[i, 1]
  target_chr  <- sig_snps[i, 2]
  target_pos  <- sig_snps[i, 3]
  trait_name <- sig_snps[i, 10]

  region_hits <- lead_coloc_snps %>% filter(Target_SNP == rsid & GWAS_Trait == trait_name)

  if(nrow(region_hits) == 0) {
    message(paste0("\n[", i, "/", nrow(sig_snps), "] Skipping ", rsid, " (No hits > 0.75)"))
    next
  }
  
  trait_gcst <- gwas_info$name[gwas_info$Disease == trait_name]
  file_gwas <- paste0("Data_processed/", rsid, "_", trait_gcst, "_", p_threshold, ".rds")
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
    
    file_eqtl <- paste0("Data_processed/eqtl_", rsid, "_", gene_name, ".rds")
    
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
                            trait_name, shared_positions, target_pos, rsid, "red", window_min, window_max)
  
    p2 <- manhattan_plot_fx(filter(df_plot, Dataset == paste0(gene_name, " eQTL")), 
                            paste0(gene_name, " eQTL"), shared_positions, target_pos, rsid, "blue", window_min, window_max)

    combined_plot <- (p1 / p2) + 
      plot_annotation(
        title = paste0("Manhattan Plot: ", rsid, " (Chromosome ", target_chr, ")"),
        subtitle = paste0("PP_H4: ", pp4)
      )
    
    safe_rsid <- gsub(":", "_", rsid)
    plot_filename <- file.path(out_dir, paste0("Plot_", safe_rsid, "_", trait_name, "_vs_", gene_name, "_eQTL.png"))
  
    ggsave(plot_filename, plot = combined_plot, width = 10, height = 8, dpi = 300)
    message(paste("  -> Saved:", plot_filename))
  }
}

message("\nAll plotting complete.")