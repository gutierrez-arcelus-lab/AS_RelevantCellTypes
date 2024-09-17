#!/bin/bash
#SBATCH --job-name=ldsc_seg
#SBATCH --output=out_ldsc_seg.txt
#SBATCH --error=err_ldsc_seg.txt
#SBATCH --partition=bch-compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Marcos.ChinasHernandez@childrens.harvard.edu
#SBATCH --time=01:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=4 

source /programs/biogrids.shrc
module load anaconda2
sh balanced_groups_tstat/ldsc_commands_balanced_groups_tstat.sh
sh balanced_groups_tstat_correctionbystim/ldsc_commands_balanced_groups_tstat_correctionbystim.sh