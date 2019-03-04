#!/bin/bash
#
#SBATCH --job-name=Zymbo_analysis
#SBATCH --output=Zymbo_analysis.out
#SBATCH --error=Zymbo_analysis.err
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --mem=30G


module load R/3.5.1

R -e "rmarkdown::render('Zymbo_analysis.Rmd')"
