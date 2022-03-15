#!/bin/bash
#SBATCH --job-name=getrs_info
#SBATCH --output=getrs_%A_%a.out
#SBATCH --error=getrs_%A_%a.err
#SBATCH --partition=bch-compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Marcos.ChinasHernandez@childrens.harvard.edu
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=8 

source /programs/biogrids.shrc

ldsc.py 


/lab-share/IM-Gutierrez-e2/Public/LDSC-Calderon/Calderon_baseline_45W
