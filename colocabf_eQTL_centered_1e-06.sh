#!/bin/bash
#SBATCH --job-name=coloceqtl_recentering
#SBATCH --account=mrc-bsu2-sl2-cpu
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/coloceqtl_recentering_%j.out
#SBATCH --error=logs/coloceqtl_recentering_%j.err
#SBATCH --partition=icelake
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xz418@cam.ac.uk

module load R/4.5
Rscript colocabf_eQTL_centered_1e-06.R
