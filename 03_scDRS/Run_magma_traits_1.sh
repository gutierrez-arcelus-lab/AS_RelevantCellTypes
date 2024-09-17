#!/bin/bash
#SBATCH --job-name=run_magma_traits_1
#SBATCH --output=run_magma_traits_1_out_%a.txt
#SBATCH --error=run_magma_traits_1_err_%a.txt
#SBATCH --partition=bch-compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Marcos.ChinasHernandez@childrens.harvard.edu
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=2
#SBATCH --array=1-18

source /programs/biogrids.shrc

commands=(
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/AS_ImmunoChip.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/AS_ImmunoChip'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/AS_finngen.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/AS_finngen'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/AS_finngen_strict.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/AS_finngen_strict'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/AS_immunochip_removed_MHC_24MB_to_35MB.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/AS_immunochip_removed_MHC_24MB_to_35MB'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/AS_panUK.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/AS_panUK'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/AdultOnset.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/AdultOnset'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/AllergyEczema.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/AllergyEczema'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/Alzheimer.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/Alzheimer'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/Arthritis.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/Arthritis'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/ChildOnset.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/ChildOnset'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/Crohns_Disease.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/Crohns_Disease'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/Height.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/Height'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/IBD.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/IBD'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/Lupus.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/Lupus'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/New_Lupus.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/New_Lupus'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/PBC.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/PBC'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/UKBAsthma.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/UKBAsthma'
  './magma_v1.10/magma --bfile ./magma_v1.10/g1000_eur/g1000_eur --pval ./magma_v1.10/Ulcerative_Colitis.pval use='SNP,P' ncol='N' --gene-annot ./magma_v1.10/out/step1.genes.annot --out ./magma_v1.10/out/step2/Ulcerative_Colitis'
)

cmd=${commands[SLURM_ARRAY_TASK_ID-1]}
echo "Running command: $cmd"
eval $cmd

