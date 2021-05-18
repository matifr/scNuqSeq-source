# Single nuq RNA-seq documentation

## File structure

1. <mark>My_initial_analysis_high_var.ipynb</mark>

2. Cell_Identity_analysis_high_var.ipynb

3. Microglia_subcluster.ipynb
   
   1. Microglia_subcluster_BS_.ipynb

4. Astrocyte_subcluster.ipynb

5. oligo_subcluster.ipynb

6. runMAST_clusters.R

7. runMAST.R

8. runMAST_condition.R

9. 

## 1. My_initial_analysis_high_var

This is the file where the preprocessing and filtering of dead cells is happening. 

**Output**: `adata_final_nov.h5ad`

Notes: Need to add doublte filtering here

## 2. Cell_Identity_analysis_high_var

This is the script where I identify the different broad clusters based on base markers. 

1. Immune cells/Microglia

2. Glia

3. Neurons

4. Epithelial cells

5. Vascular/Fibro cells

**Output**: `adata_final_high_var_broad_types`



For each broad cluster, I subcluster into its cell types    

- <mark>Immune cells/Microglia:</mark> **Output**: `adata_final_high_var_immune_cells`

- <mark>Glia:</mark> **Output**: `adata_final_high_var_glia_cells`

- <mark>Neurons</mark> **Output**: `adata_final_high_var_neuron_cells`



Final cluster that has all the information in one AnnData

`adata_final_high_var_final_clusters_updated`



## 3. Microglia_subcluster

- File that inputs `adata_final_high_var_final_clusters_updated`

- Subsets the microlgia, does all the preprocessing, nornalization and clustering.

- Identifies cluster <mark>F131A</mark> **Output**:`adata_high_var_Microglia_clusters`

- Excludes F131A and reclusters the rest of microlgia and **Output**:`adata_high_var_Microglia_ONLY_clusters`



### 3a. Microglia_subcluster_BS

**Input**:`adata_high_var_Microglia_ONLY_clusters`

This script subclusters the <mark>brainstem</mark> microlgia, does all the preprocessing, nornalization and clustering. Conducts DEA, heatmaps for the specific microglia subsets i.e. APOE high/low

- **Output**:`adata_high_var_Microglia_Brainstem_clusters`




















