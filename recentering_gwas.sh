#!/bin/bash
#SBATCH --job-name=gwas_recentering
#SBATCH --account=mrc-bsu2-sl2-cpu
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/recentering_gwas_%j.out
#SBATCH --error=logs/recentering_gwas_%j.err
#SBATCH --partition=icelake
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xz418@cam.ac.uk

module load R/4.5
Rscript recentering_gwas.R
