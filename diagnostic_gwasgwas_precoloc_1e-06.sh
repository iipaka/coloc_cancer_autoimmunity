#!/bin/bash
#SBATCH --job-name=plot_diagnostic_precoloc
#SBATCH --account=mrc-bsu2-sl2-cpu
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/diagnostic_gwasgwas_precoloc_1e-06_%j.out
#SBATCH --error=logs/diagnostic_gwasgwas_precoloc_1e-06_%j.err
#SBATCH --partition=icelake
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xz418@cam.ac.uk

module load R/4.3.1-gcc
Rscript diagnosticplot_gwasgwas_precoloc_1e-06.R
