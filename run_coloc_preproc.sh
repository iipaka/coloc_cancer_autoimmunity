#!/bin/bash
#SBATCH --job-name=coloc_job
#SBATCH --account=mrc-bsu2-sl2-cpu
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:45:00
#SBATCH --mem=32G
#SBATCH --output=logs/coloc_preproc_log_%j.out
#SBATCH --error=logs/coloc_preproc_error_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xz418@cam.ac.uk


module load R/4.5
Rscript clean_gwas.R
