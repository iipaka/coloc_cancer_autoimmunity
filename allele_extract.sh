#!/bin/bash
#SBATCH --job-name=allele_extract
#SBATCH --account=mrc-bsu2-sl2-cpu
#SBATCH --time=00:10:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/allele_extract_%j.out
#SBATCH --error=logs/allele_extract_%j.err
#SBATCH --partition=icelake
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xz418@cam.ac.uk

module load R/4.3.1-icelake
Rscript extract_alleles.R
