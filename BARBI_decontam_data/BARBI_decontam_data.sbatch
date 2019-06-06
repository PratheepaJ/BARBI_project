#!/bin/bash
#
#SBATCH --job-name=BARBI_decontam_data
#SBATCH --output=BARBI_decontam_data.out
#SBATCH --error=BARBI_decontam_data.err
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --mem=30G


module load R/3.5.1

R -e "rmarkdown::render('BARBI_decontam_data.Rmd')"
