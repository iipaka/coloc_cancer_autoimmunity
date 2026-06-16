#!/bin/bash
#SBATCH --job-name=eQTL_rds_extraction
#SBATCH --account=mrc-bsu2-sl2-cpu
#SBATCH --time=00:10:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/eQTL_rds_extraction_%j.out
#SBATCH --error=logs/eQTL_rds_extraction_%j.err
#SBATCH --partition=icelake
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xz418@cam.ac.uk

module load R/4.5
module load htslib/1.14/gcc/wwc5wqv5
Rscript eQTL_rds_extraction_for_coloc.R
