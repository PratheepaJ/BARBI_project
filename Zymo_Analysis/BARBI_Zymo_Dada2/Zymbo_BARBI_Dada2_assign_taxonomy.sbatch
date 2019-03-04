#!/bin/bash
#
#SBATCH --job-name=Zymbo_BARBI_Dada2_assign_taxonomy
#SBATCH --output=Zymbo_BARBI_Dada2_assign_taxonomy.out
#SBATCH --error=Zymbo_BARBI_Dada2_assign_taxonomy.err
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --mem=30G


module load R/3.5.1

R -e "rmarkdown::render('Zymbo_BARBI_Dada2_assign_taxonomy.Rmd')"
