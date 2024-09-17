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


cd /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim
sh /lab-share/IM-Gutierrez-e2/Public/Dany/PYTHON-BASH_SCRIPTS/window-LD2.sh -i /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/BCells_mixedmodel.bed -w 225
sh /lab-share/IM-Gutierrez-e2/Public/Dany/PYTHON-BASH_SCRIPTS/window-LD2.sh -i /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/DC1_mixedmodel.bed -w 225
sh /lab-share/IM-Gutierrez-e2/Public/Dany/PYTHON-BASH_SCRIPTS/window-LD2.sh -i /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/DC2_mixedmodel.bed -w 225
sh /lab-share/IM-Gutierrez-e2/Public/Dany/PYTHON-BASH_SCRIPTS/window-LD2.sh -i /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/Monocytes_mixedmodel.bed -w 225
sh /lab-share/IM-Gutierrez-e2/Public/Dany/PYTHON-BASH_SCRIPTS/window-LD2.sh -i /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/NKs_mixedmodel.bed -w 225
sh /lab-share/IM-Gutierrez-e2/Public/Dany/PYTHON-BASH_SCRIPTS/window-LD2.sh -i /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/Plasma_mixedmodel.bed -w 225
sh /lab-share/IM-Gutierrez-e2/Public/Dany/PYTHON-BASH_SCRIPTS/window-LD2.sh -i /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/TCells_mixedmodel.bed -w 225
echo BCells_mixedmodel /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/225/BCells_mixedmodel_225w.,/lab-share/IM-Gutierrez-e2/Public/References/LDSC-Calderon/Control_w225/Control_w225. > /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/balanced_groups_tstat_correctionbystim.ldcts
echo DC1_mixedmodel /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/225/DC1_mixedmodel_225w.,/lab-share/IM-Gutierrez-e2/Public/References/LDSC-Calderon/Control_w225/Control_w225. >> /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/balanced_groups_tstat_correctionbystim.ldcts
echo DC2_mixedmodel /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/225/DC2_mixedmodel_225w.,/lab-share/IM-Gutierrez-e2/Public/References/LDSC-Calderon/Control_w225/Control_w225. >> /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/balanced_groups_tstat_correctionbystim.ldcts
echo Monocytes_mixedmodel /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/225/Monocytes_mixedmodel_225w.,/lab-share/IM-Gutierrez-e2/Public/References/LDSC-Calderon/Control_w225/Control_w225. >> /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/balanced_groups_tstat_correctionbystim.ldcts
echo NKs_mixedmodel /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/225/NKs_mixedmodel_225w.,/lab-share/IM-Gutierrez-e2/Public/References/LDSC-Calderon/Control_w225/Control_w225. >> /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/balanced_groups_tstat_correctionbystim.ldcts
echo Plasma_mixedmodel /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/225/Plasma_mixedmodel_225w.,/lab-share/IM-Gutierrez-e2/Public/References/LDSC-Calderon/Control_w225/Control_w225. >> /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/balanced_groups_tstat_correctionbystim.ldcts
echo TCells_mixedmodel /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/225/TCells_mixedmodel_225w.,/lab-share/IM-Gutierrez-e2/Public/References/LDSC-Calderon/Control_w225/Control_w225. >> /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/01_LDSC-SEG/LDSC_seg/balanced_groups_tstat_correctionbystim/balanced_groups_tstat_correctionbystim.ldcts
