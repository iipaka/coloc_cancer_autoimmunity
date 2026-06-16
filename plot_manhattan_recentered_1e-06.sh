#!/bin/bash
#SBATCH --job-name=plot_manhattan_eQTL_recentered
#SBATCH --account=mrc-bsu2-sl2-cpu
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/plot_manhattan_eQTL_recentered_1e-06_%j.out
#SBATCH --error=logs/plot_manhattan_eQTL_recentered_1e-06_%j.err
#SBATCH --partition=icelake
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xz418@cam.ac.uk

module load R/4.5
Rscript Manhattan_eQTL_centered_1e-06.R
