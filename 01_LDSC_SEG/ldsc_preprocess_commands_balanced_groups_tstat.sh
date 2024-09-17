#!/bin/bash
#SBATCH --job-name=ldsc_preprocess
#SBATCH --output=out_ldsc_preprocess.txt
#SBATCH --error=err_ldsc_preprocess.txt
#SBATCH --partition=bch-compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Marcos.ChinasHernandez@childrens.harvard.edu
#SBATCH --time=01:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=4 

source /programs/biogrids.shrc
module load anaconda2


cd /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LDSC_seg/balanced_groups_tstat
sh /lab-share/IM-Gutierrez-e2/Public/Dany/PYTHON-BASH_SCRIPTS/window-LD2.sh -i /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LDSC_seg/balanced_groups_tstat/ -w 225
echo NA /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LDSC_seg/balanced_groups_tstat/225/NA_225w.,/lab-share/IM-Gutierrez-e2/Public/References/LDSC-Calderon/Control_w225/Control_w225. > /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LDSC_seg/balanced_groups_tstat/balanced_groups_tstat.ldcts
echo  /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LDSC_seg/balanced_groups_tstat/225/_225w.,/lab-share/IM-Gutierrez-e2/Public/References/LDSC-Calderon/Control_w225/Control_w225. >> /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LDSC_seg/balanced_groups_tstat/balanced_groups_tstat.ldcts
