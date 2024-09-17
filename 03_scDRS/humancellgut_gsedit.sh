#!/bin/bash
#SBATCH --job-name=immunecellgut_traits
#SBATCH --output=out/output_immunecellgut_traits_gsedit_%a.txt
#SBATCH --error=out/error-immunecellgut_traits_%a.txt
#SBATCH --partition=bch-compute-pe
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Marcos.ChinasHernandez@childrens.harvard.edu
#SBATCH --time=07:30:00
#SBATCH --mem=280G
#SBATCH --ntasks=3
#SBATCH --array=1-2%1
source /programs/biogrids.shrc

scdrs_dir="/lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis"
project="humancellgut_traits_gsedit"
cov="/lab-share/IM-Gutierrez-e2/Public/References/humancellatlas/humancellgut_run_cov.cov"
my_h5ad="/lab-share/IM-Gutierrez-e2/Public/References/humancellatlas/Full_obj_raw_counts_nosoupx_v2.h5ad"
adj_by="Integrated_05"
gs_file_dir="gs_file_gutcellatlas"

#commands=(AS_ImmunoChip	AS_finngen	AS_finngen_strict	AS_immunochip_removed_MHC_24MB_to_35MB	AS_panUK	Alzheimer	Arthritis	Height	IBD	UKBAsthma	Crohns_Disease	Lupus	New_Lupus	PBC	Ulcerative_Colitis)
commands=(AdultOnset	ChildOnset	AS_finngen_strict	AS_immunochip_removed_MHC_24MB_to_35MB	AS_panUK	Alzheimer	Arthritis	Height	IBD	UKBAsthma	Crohns_Disease	Lupus	New_Lupus	PBC	Ulcerative_Colitis)


trait=${commands[SLURM_ARRAY_TASK_ID-1]}
#echo $trait >> begin_immunecellgut_traits.txt

source "/lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_scripts/functions_scdrs.sh"
create_dirs
scdrs_cov
add_scdrs_meta "cov"
### Run on a next cycle and comment scdrs_cov and add_scdrs_meta ### 
#scdrs_downstream "cov" 


#scdrs_default
#scdrs_adjprop
#scdrs_adjprop_cov
#scdrs_filtfalse
#add_scdrs_meta "default"


