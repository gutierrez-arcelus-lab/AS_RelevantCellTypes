#!/bin/bash
#SBATCH --job-name=run_magma_traits_3
#SBATCH --output=run_magma_traits_3_out_%a.txt
#SBATCH --error=run_magma_traits_3_err_%a.txt
#SBATCH --partition=bch-compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Marcos.ChinasHernandez@childrens.harvard.edu
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1

source /programs/biogrids.shrc
module load singularity

singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/AS_ImmunoChip.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/AS_ImmunoChip_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/AS_finngen.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/AS_finngen_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/AS_finngen_strict.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/AS_finngen_strict_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/AS_immunochip_removed_MHC_24MB_to_35MB.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/AS_immunochip_removed_MHC_24MB_to_35MB_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/AS_panUK.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/AS_panUK_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/AdultOnset.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/AdultOnset_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/AllergyEczema.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/AllergyEczema_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/Alzheimer.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/Alzheimer_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/Arthritis.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/Arthritis_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/ChildOnset.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/ChildOnset_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/Crohns_Disease.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/Crohns_Disease_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/Height.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/Height_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/IBD.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/IBD_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/Lupus.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/Lupus_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/New_Lupus.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/New_Lupus_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/PBC.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/PBC_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/UKBAsthma.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/UKBAsthma_scDRS_input_magma.tsv --weight zscore --n-max 1000
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs munge-gs --out-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis/gs_file_grch37/Ulcerative_Colitis.gs --zscore-file /lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/magma_v1.10/gs_file/Ulcerative_Colitis_scDRS_input_magma.tsv --weight zscore --n-max 1000


