# ==============================================================================
# Purpose: Generate Manhattan plots (GWAS <> eQTL) on Centred Data (PRE-COLOC)
# ==============================================================================
library(ggplot2)
library(dplyr)
library(data.table)
library(patchwork)

p_threshold <- "1e-06"

out_dir <- paste0("Output/Plots_eQTL_", p_threshold, "_centered")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sig_snps <- read.csv("Output/recentered_targets.csv")
gwas_info <- read.csv("Data_raw/gwas_summarystats.csv")

# ==============================================================================
# PLOTTING FUNCTION
# ==============================================================================

manhattan_theme <- theme_bw(base_family = "Arial") + 
  theme(
    text = element_text(family = "Arial"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 10, face = "plain")
  )

gwas_threshold <- -log10(as.numeric(p_threshold))

manhattan_plot_fx <- function(data, title, shared_pos, target_pos, dot_colour, x_min, x_max) {
  
  ggplot(data, aes(x = POS / 1e6, y = logP)) +
    #Vertical line anchoring the centre
    geom_vline(xintercept = target_pos / 1e6, linetype = "dashed", alpha = 0.5) +
    
    #Shared points
    geom_point(data = filter(data, POS %in% shared_pos), 
               aes(color = "Common to both"), alpha = 0.7, size = 1.75, stroke = 0.5) +
    
    #Base points
    geom_point(data = filter(data, !(POS %in% shared_pos)), 
               aes(color = Dataset), alpha = 0.6, size = 1.75, stroke = 0.5) +
    
    #Scale and labels
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
  
  clean_snp <- sub(":.*", "", sig_snps$New_Centre_SNP[i])
  target_snp <- sig_snps$New_Centre_SNP[i]
  target_chr <- sig_snps$CHR[i]
  target_pos <- sig_snps$New_Centre_POS[i]
  trait_name <- sig_snps$Stronger_Trait[i]
  
  trait_gcst <- gwas_info$name[gwas_info$Disease == trait_name]
  file_gwas <- paste0("Data_processed/", target_snp, "_", trait_gcst, "_", p_threshold, "_centered.rds")
  
  data_gwas <- readRDS(file_gwas)
  df_gwas <- data.frame(
    POS = data_gwas$position,
    P = data_gwas$pvalues,
    Dataset = trait_name 
  )
  
  eqtl_pattern <- paste0("^eqtl_", target_snp, "_.*_centered\\.rds$")
  eqtl_files <- list.files("Data_processed", pattern = eqtl_pattern, full.names = TRUE)
  
  for(eqtl_f in eqtl_files) {
    
    file_base <- gsub("_centered\\.rds$", "", basename(eqtl_f))
    gene_name <- strsplit(file_base, "_")[[1]][3]
    
    message(paste0("  -> Plotting hit: ", trait_name, " vs ", gene_name, " eQTL"))
    
    data_eqtl <- readRDS(eqtl_f)
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
        title = paste0("Manhattan Plot (centered): ", clean_snp, " (Chromosome ", target_chr, ")"),
        theme = theme(plot.title = element_text (size = 22, family = "Arial")))
    
    safe_gene <- gsub(":", "-", gene_name)
    plot_filename <- file.path(out_dir, paste0("Plot_", clean_snp, "_", trait_name, "_vs_", safe_gene, "_eQTL.png"))
    
    ggsave(plot_filename, plot = combined_plot, width = 10, height = 8, dpi = 300)
  }
}

message("\nAll plotting complete.")