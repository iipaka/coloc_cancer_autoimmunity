# ==============================================================================
# Purpose: Generate regional Manhattan plots for target SNPs
# ==============================================================================
library(ggplot2)
library(dplyr)
library(data.table)
library(patchwork)

p_threshold <- 1e-6

out_dir <- paste0("Output/Plots_", p_threshold)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sig_snps <- read.csv(paste0("Output/sig_snps_threshold_", p_threshold, ".csv"))
gwas_info <- read.csv("Data_raw/gwas_summarystats.csv")

# ==============================================================================
# PLOTTING
# ==============================================================================

manhattan_theme <- theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 10, face = "plain")
  )

gws_threshold <- -log10(p_threshold)

manhattan_plot_fx <- function(data, title, shared_pos, target_pos, target_name, dot_colour, x_min, x_max) {
  
  ggplot(data, aes(x = POS / 1e6, y = logP)) +
    
    #Vertical line for target
    geom_vline(xintercept = target_pos / 1e6, linetype = "dashed", alpha = 0.5) +
    
    #Shared points
    geom_point(data = filter(data, POS %in% shared_pos), 
               aes(color = "Common to both"), alpha = 0.7, size = 1.75, stroke = 0.5) +
    
    #Base points
    geom_point(data = filter(data, !(POS %in% shared_pos)), 
               aes(color = Disease), alpha = 0.6, size = 1.75, stroke = 0.5) +
    
    #Target SNP highlighting
    geom_point(data = filter(data, POS == target_pos), 
               shape = 1, size = 4, color = "black", stroke = 1) +
    
    #Target SNP text
    geom_text(data = filter(data, POS == target_pos), 
              aes(label = target_name), 
              hjust = -0.2, vjust = 0.5, fontface = "bold", size = 3.5) +
    
    #Significance line and text
    geom_hline(yintercept = gws_threshold, linetype = "dotted", color = "black") +
    annotate("text", x = Inf, y = gws_threshold, label = "p-value threshold", 
             vjust = -0.5, hjust = 1.05, color = "black", size = 3.5) +
    
    #Scale and Labels
    scale_color_manual(
      name = "SNP", 
      breaks = c(title, "Common to both"), 
      values = setNames(c("pink", dot_colour), c("Common to both", title))
    ) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    labs(title = title, x = expression("BP38 / 10"^6), y = expression("-log"[10]*"(P-value)")) +
    manhattan_theme +
    
    theme(
      plot.title = element_text(size = 18),
      axis.title.x = element_text(size = 16, margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, margin = margin(r = 10)),
      axis.text.x = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
}

# ==============================================================================
# LOOP THROUGH SNPs
# ==============================================================================

for(i in 1:nrow(sig_snps)) {
  
  #Extract relevant columns for mapping
  rsid <- sig_snps[i, 1]
  target_chr <- sig_snps[i, 2]
  target_pos <- sig_snps[i, 3]
  cancer_pval <- formatC(sig_snps[i, 9], format = "e", digits = 2) 
  immune_name <- sig_snps[i, 10]
  immune_pval <- formatC(sig_snps[i, 13], format = "e", digits = 2)
  cancer_name <- sig_snps[i, 16]
  
  message(paste0("\n[", i, "/", nrow(sig_snps), "] Plotting ", rsid, " (", cancer_name, " vs ", immune_name, ")"))
  
  #Find the correct GCST IDs 
  cancer_gcst <- gwas_info$name[gwas_info$Disease == cancer_name]
  immune_gcst <- gwas_info$name[gwas_info$Disease == immune_name]
  
  #Construct file paths
  file_cancer <- paste0("Data_processed/", rsid, "_", cancer_gcst, "_", p_threshold, ".rds")
  file_immune <- paste0("Data_processed/", rsid, "_", immune_gcst, "_", p_threshold, ".rds")
  
  #Load data
  data_cancer <- readRDS(file_cancer)
  data_immune <- readRDS(file_immune)
  
  df_cancer <- data.frame(
    POS = data_cancer$position,
    P = data_cancer$pvalues,
    Disease = cancer_name 
  )
  
  df_immune <- data.frame(
    POS = data_immune$position,
    P = data_immune$pvalues,
    Disease = immune_name 
  )
  
  #Data processing for plot
  shared_positions <- intersect(df_cancer$POS, df_immune$POS)
  
  df_plot <- rbind(df_cancer, df_immune)
  df_plot$logP <- -log10(df_plot$P + 1e-300) 
  
  window_min <- min(df_plot$POS, na.rm = TRUE) / 1e6
  window_max <- max(df_plot$POS, na.rm = TRUE) / 1e6
  
  #Generate plots
  p1 <- manhattan_plot_fx(filter(df_plot, Disease == cancer_name), 
                            cancer_name, shared_positions, target_pos, rsid, "red", window_min, window_max)
  
  p2 <- manhattan_plot_fx(filter(df_plot, Disease == immune_name), 
                            immune_name, shared_positions, target_pos, rsid, "blue", window_min, window_max)
  
  #Combine 
  combined_plot <- (p1 / p2) + 
    plot_annotation(title = paste0("Manhattan Plot: ", rsid, " (Chromosome ", target_chr, ")"),
                    theme = theme(plot.title = element_text (size = 22)))
  
  #Save
  safe_rsid <- gsub(":", "_", rsid)
  plot_filename <- file.path(out_dir, paste0("Plot_", safe_rsid, "_", cancer_name, "_vs_", immune_name, ".png"))
  
  ggsave(plot_filename, plot = combined_plot, width = 10, height = 8, dpi = 300)
  message(paste("  -> Saved:", plot_filename))
}

message("\nAll plotting complete.")