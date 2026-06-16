#!/bin/bash
#SBATCH --job-name=coloc_15_regions
#SBATCH --account=mrc-bsu2-sl2-cpu
#SBATCH --partition=icelake
#SBATCH --output=logs/colocabf_1e-06_%j.out
#SBATCH --error=logs/colocabf_1e-06_%j.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xz418@cam.ac.uk

module load R/4.5
Rscript colocabf_all_1e-06.R

