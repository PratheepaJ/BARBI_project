#!/bin/bash
#
#SBATCH --job-name=HealthyFlossing_BARBI
#SBATCH --output=HealthyFlossing_BARBI.out
#SBATCH --error=HealthyFlossing_BARBI.err
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mem=20G

module load R/3.5.1

R -e "rmarkdown::render('HealthyFlossing_BARBI.Rmd')"
