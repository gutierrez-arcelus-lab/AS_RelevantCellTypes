#!/bin/bash
#SBATCH --job-name=LMM_Balanced
#SBATCH --partition=bch-compute
#SBATCH -n 10
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Marcos.ChinasHernandez@childrens.harvard.edu
#SBATCH --output=/lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LMM_LDSC_seg/LMM_Balanced_%A.out
#SBATCH --error=/lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LMM_LDSC_seg/LMM_Balanced_%A.err
#SBATCH --time=10:00:00


source /programs/biogrids.shrc
cd /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LMM_LDSC_seg
export R_X=3.6.2
# This script uses a modified version of 
Rscript /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LMM_LDSC_seg/peak_selection_balanced_groups.R \
-o "/lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LMM_LDSC_seg/balanced_groups" \
-p "/lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LMM_LDSC_seg/balanced_groups.csv" 

