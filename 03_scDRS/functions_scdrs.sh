#!/bin/bash
create_dirs () {
if [  -d "${scdrs_dir}/${project}" ]; then echo "Files exits"; return 1; fi
mkdir  "${scdrs_dir}/${project}"
mkdir  "${scdrs_dir}/${project}/scdrs_results/"

mkdir  "${scdrs_dir}/${project}/scdrs_results/default"
mkdir  "${scdrs_dir}/${project}/scdrs_results/default/Downstream"

mkdir  "${scdrs_dir}/${project}/scdrs_results/cov"
mkdir  "${scdrs_dir}/${project}/scdrs_results/cov/Downstream"

mkdir  "${scdrs_dir}/${project}/scdrs_results/adjprop"
mkdir  "${scdrs_dir}/${project}/scdrs_results/adjprop/Downstream"

mkdir  "${scdrs_dir}/${project}/scdrs_results/adjprop_cov"
mkdir  "${scdrs_dir}/${project}/scdrs_results/adjprop_cov/Downstream"

mkdir  "${scdrs_dir}/${project}/scdrs_results/cov_filtfalse"
mkdir  "${scdrs_dir}/${project}/scdrs_results/cov_filtfalse/Downstream"

mkdir  "${scdrs_dir}/${project}/scdrs_results/filtfalse"
mkdir  "${scdrs_dir}/${project}/scdrs_results/filtfalse/Downstream"

}
export SINGULARITY_CACHEDIR='/temp_work/ch229505/'
export SINGULARITY_DOCKER_USERNAME=marcoschinas
export SINGULARITY_DOCKER_PASSWORD='smQ)rKL?~34aEBJ'

scdrs_default () {
#module load singularity
cd "/lab-share/IM-Gutierrez-e2/Public/" 
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs compute-score \
--h5ad_file ${my_h5ad} \
--h5ad_species human \
--gs_file "${scdrs_dir}/${gs_file_dir}/${trait}.gs" \
--gs_species human \
--out_folder ${scdrs_dir}/${project}/scdrs_results/default 
}

#Requires additionally covariate file
scdrs_mm_default () {
#module load singularity
cd "/lab-share/IM-Gutierrez-e2/Public/" 
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs compute-score \
--h5ad_file ${my_h5ad} \
--h5ad_species mouse \
--gs_file "${scdrs_dir}/${gs_file_dir}/${trait}.gs" \
--gs_species mouse \
--out_folder ${scdrs_dir}/${project}/scdrs_results/default
}
#Requires additionally covariate file

scdrs_cov () {
#module load singularity
cd "/lab-share/IM-Gutierrez-e2/Public/" 
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs compute-score \
--h5ad_file ${my_h5ad} \
--h5ad_species human \
--cov_file ${cov} \
--gs_file "${scdrs_dir}/${gs_file_dir}/${trait}.gs" \
--gs_species human \
--out_folder ${scdrs_dir}/${project}/scdrs_results/cov
}

scdrs_mm_cov () {
#module load singularity
cd "/lab-share/IM-Gutierrez-e2/Public/" 
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs compute-score \
--h5ad_file ${my_h5ad} \
--h5ad_species mouse \
--cov_file ${cov} \
--gs_file "${scdrs_dir}/${gs_file_dir}/${trait}.gs" \
--gs_species mouse \
--out_folder ${scdrs_dir}/${project}/scdrs_results/cov
}

scdrs_cov_5k () {
#module load singularity
cd "/lab-share/IM-Gutierrez-e2/Public/" 
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs compute-score \
--h5ad_file ${my_h5ad} \
--h5ad_species human \
--cov_file ${cov} \
--gs_file "${scdrs_dir}/${gs_file_dir}/${trait}.gs" \
--gs_species human \
--out_folder ${scdrs_dir}/${project}/scdrs_results/cov \
--n-ctrl 5000
}

scdrs_adjprop_cov () {
#module load singularity
cd "/lab-share/IM-Gutierrez-e2/Public/" 

singularity run docker://marcoschinas/scdrs:1.0.2 scdrs compute-score \
--h5ad_file ${my_h5ad} \
--h5ad_species human \
--cov_file ${cov} \
--gs_file "${scdrs_dir}/${gs_file_dir}/${trait}.gs" \
--gs_species human \
--adj_prop ${adj_by} \
--out_folder ${scdrs_dir}/${project}/scdrs_results/adjprop_cov
}

### Requires column to adjust proportions ###
scdrs_adjprop () {
#module load singularity
cd "/lab-share/IM-Gutierrez-e2/Public/" 
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs compute-score \
--h5ad_file ${my_h5ad} \
--h5ad_species human \
--gs_file "${scdrs_dir}/${gs_file_dir}/${trait}.gs" \
--gs_species human \
--adj_prop ${adj_by} \
--out_folder ${scdrs_dir}/${project}/scdrs_results/adjprop/ 
}
# Perform scdrs analysis withouth filtering
scdrs_filtfalse () { 
#module load singularity 
cd "/lab-share/IM-Gutierrez-e2/Public/" 
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs compute-score \
--h5ad_file ${my_h5ad} \
--h5ad_species human \
--gs_file "${scdrs_dir}/${gs_file_dir}/${trait}.gs" \
--gs_species human \
--out_folder ${scdrs_dir}/${project}/scdrs_results/filtfalse \
--flag_filter_data FALSE
}

### Perform scdrs downstream analysis, specify folder to perform the analysis, otherwise will take the full scores from default
scdrs_downstream () { 
folder=${1:-default}
#module load singularity 
cd "/lab-share/IM-Gutierrez-e2/Public/" 
singularity run docker://marcoschinas/scdrs:1.0.2 scdrs perform-downstream \
--h5ad-file ${my_h5ad} \
--score-file "${scdrs_dir}/${project}/scdrs_results/${folder}/${trait}.full_score.gz" \
--out-folder ${scdrs_dir}/${project}/scdrs_results/${folder}/Downstream \
--flag-filter-data True \
--gene-analysis \
--group_analysis ${adj_by} \
--flag-raw-count True
}

# Add scdrs results using # 
add_res_scdrs () {
(
source /programs/biogrids.shrc
export R_X=4.1
replace=".h5ad"
replacewith=".h5Seurat"
my_h5Seurat=$(sed 's/'$replace'$/'$replacewith'/g' <<<"$my_h5ad")
replacewith="_metadata.txt"
meta=$(sed 's/'$replace'$/'$replacewith'/g' <<<"$my_h5ad")
folder=${1:-default}
split=${2:-FALSE}

Rscript '/lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_scripts/add_res_scdrs.R' \
-r ${scdrs_dir}/${project}/scdrs_results/$folder \
-i  ${my_h5Seurat} -c ${adj_by} -o "${scdrs_dir}/${project}/scdrs_results/$folder/Summary" -s ${split}
)
}
# Add scdrs using metadata #
add_scdrs_meta () {
(
source /programs/biogrids.shrc
export R_X=4.1
replace=".h5ad"
replacewith=".h5Seurat"
my_h5Seurat=$(sed 's/'$replace'$/'$replacewith'/g' <<<"$my_h5ad")
replacewith="_metadata.txt"
meta=$(sed 's/'$replace'$/'$replacewith'/g' <<<"$my_h5ad")

echo $meta
folder=${1:-default}
split=${2:-FALSE}
Rscript '/lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_scripts/add_scdrs_to_meta.R' \
-r ${scdrs_dir}/${project}/scdrs_results/$folder \
-i  ${meta} \
-c ${adj_by} -o "${scdrs_dir}/${project}/scdrs_results/$folder/Summary" -s ${split}
)
}

