# Repository for running colocalisation analysis between cancer and immune-mediated disease

This repository contains the R scripts and analytical pipeline used in the dissertation: *"Unravelling Opposing Effects of Genetic Variants in Cancer and Autoimmunity to Inform Immunotherapy"*. The study investigates the shared genetic architecture between cancer and immune-mediated diseases (IMDs) using statistical colocalisation (GWAS-GWAS and GWAS-eQTL). R version 4.3.1 was used to run the analyses for the dissertation.

## Directory Structure

To replicate this analysis, the local workspace should match the following directory structure before running the scripts:

```text
├── Data_raw/            # Place downloaded GWAS summary statistics (.gz), the Chen et al. (2025) supplementary table, and the gwas_summarystats.csv mapping file here.
├── Data_processed/      # Output directory for intermediate formatted files and .rds objects required by the coloc package.
├── Output/              # Output directory for final colocalisation results, logs, and generated plots.
```

### Central Mapping File
*   `gwas_summarystats.csv`: This file is placed in `Data_raw/`. It acts as the central mapping directory for the pipeline, containing file paths to the raw `.gz` datasets, disease names, total sample size (`N`), case-to-total sample proportion (`s`), and GWAS accession numbers.

---

## Analytical Pipeline & Script Overview

The scripts are organised into folders roughly reflecting the chronological workflow of the project.

### 1. Data Collection & Formatting (`data collection/`)
*   `extract_headers.R`: Extracts and prints the column headers for each raw GWAS summary statistic file to inform downstream standardisation.
*   `BC_dataconversion.R`: Extracts and formats the Ostrogen Receptor-positive (ER+) and -negative (ER-) breast cancer summary statistics. Performs a liftover from GRCh37 to GRCh38 based on the source file headers and exports them as `.gz` files.
*   `OC_dataconversion.R`: Processes the ovarian cancer summary statistics (originally supplied as 23 separate chromosome files with an unspecified genome build). The script samples SNPs against an existing database to confirm they are GRCh37, merges the 23 files into a single dataset, performs a liftover to GRCh38, standardises the format and exports them as `.gz` files.

### 2. SNP Filtering (`snp filtering/`)
*   `filter_snps_1e-06.R`: Uses the Chen et al. (2025) supplementary data as input (`genes-3623764-supplementary.xlsx`), downloadable at https://www.mdpi.com/article/10.3390/genes16050575/s1 (also to be placed in `Data_raw/`). Filters for lead SNPs exhibiting an association strength of $p < 10^{-6}$ in both the cancer and its paired immune-mediated disease. Performs a liftover from GRCh37 to GRCh38 for the identified genomic coordinates.

### 3. Data Preparation for GWAS-GWAS Colocalisation (`data prep for gwas-gwas coloc/`)
*   `clean_gwas_1e-06.R`: Standardises summary statistics and constructs genomic windows ($\pm$ 500kb) around each filtered SNP. Uses a Unix/Linux shell command to extract specific regions from the raw GWAS files. Outputs `.rds` objects formatted to the `coloc` package's requirements. 
*   `test_coloc.R`: A diagnostic script to verify that the generated `.rds` objects are correctly formatted and ready for colocalisation.

### 4. GWAS-GWAS Colocalisation & Visualisation (`gwas-gwas coloc and plots/`)
*   `colocabf_all_1e-06.R`: Executes pairwise colocalisation on the filtered cancer/IMD pairs using the generated `.rds` objects to test for a shared causal variant.
*   `Manhattan_all_1e-06.R`: Generates regional Manhattan plot pairs to visually assess signal alignment and variant coverage across the paired traits.
*   `diagnosticplot_gwasgwas_precoloc_1e-06.R`: Generates pairwise p-value correlation plots and Z-score scatterplots to visually validate regional association patterns.
*   `diagnosticplot_gwasgwas_narrow_1e-06.R`: A post-colocalisation refinement script. For genomic regions capturing multiple distinct association signals, this script narrows the genomic window (based on visual inspection of the Manhattan plots) to exclude irrelevant secondary signals and regenerates the p-value/Z-score scatterplots.

### 5. Allele Effect Direction Verification (`allele effect direction verification/`)
*   `allele_check_x.R` *(where x = disease name)*: Investigates specific SNPs that exhibited controversial or contradictory effect directions/allele assignments in the reference meta-analysis. Checks the raw summary statistics to confirm the reference allele, effect allele, and beta coefficient signs.
*   `allele_check_cd2.R`: Performs an additional layer of validation by checking controversial Crohn's disease (CD) SNPs against an independent CD GWAS dataset.

### 6. Nearest Gene Analysis (`nearest gene analysis/`)
*   `nearest_gene_1e-06.R`: Maps the genomic loci demonstrating strong evidence of colocalisation (Posterior Probability of $H_4 > 0.75$) to their nearest protein-coding genes to postulate candidate effector targets.

### 7. Data Preparation for GWAS-eQTL Colocalisation (`data prep for gwas-eqtl coloc/`)
*   `prepare_gwas_targets_for_eqtl_coloc.R`: Identifies the statistically stronger GWAS trait (lowest p-value) within each colocalised cancer/IMD pair for eQTL colocalisation.
*   `recentering_snp.R`: Identifies the top regional GWAS SNP (strongest association) within the original 1Mb genomic window for the selected trait.
*   `recentering_gwas.R`: Re-centres the 1Mb genomic boundaries around the newly identified top regional SNP. Generates updated GWAS `.rds` files.
*   `recentering_eqtl.R`: Extracts corresponding eQTL summary statistics based on the newly updated genomic windows for all unique eGenes in the region. Generates formatted eQTL `.rds` files.

### 8. GWAS-eQTL Colocalisation & Visualisation (`gwas-eqtl coloc and plots/`)
*   `colocabf_eQTL_centered_1e-06.R`: Executes colocalisation between the recentred GWAS signals and their overlapping *cis*-eQTLs.
*   `Manhattan_eQTL_centered_1e-06.R`: Generates Manhattan plots comparing the GWAS and eQTL datasets to ensure adequate shared variant coverage.
*   `diagnosticplot_eqtlgwas_centered_precoloc_1e-06.R`: Generates p-value correlation and Z-score scatterplots for the GWAS-eQTL pairs.
*   `diagnosticplot_eqtlgwas_narrow_1e-06.R`: A post-colocalisation script that utilises a narrowed genomic window (excluding non-relevant signals) to generate refined Z-score plots, allowing for the determination of the regulatory effect direction (e.g., target gene up- or down-regulation) at that genomic window.