TrokaChatML Analysis Pipeline: From Single-Cell Data to Drug Target
Discovery
================
Michael E. Troka et al.
2025-02-25

# Introduction

TrokaChatML is an innovative computational framework for analyzing
single-cell RNA sequencing (scRNAseq) data to uncover perturbed
cell–cell communications and identify potential therapeutic targets. As
detailed in the accompanying manuscript, the pipeline integrates three
key modules:

1.  **Perturbed Cell–Cell Communication Prediction:**  
    A two-tiered noise filtration strategy is applied to remove
    artifacts inherent in scRNAseq data. First, a Gene-Wise
    Percent-Based Z Test assesses the reproducibility of gene expression
    within cell types. Then, two LightGBM machine learning models are
    executed in Python—one to filter out noise and another to determine
    which communications are truly perturbed—before this information is
    brought back into R. Notably, communication scores are computed for
    both the null distribution (derived from permutations of control
    data) and for the actual data. The scores from the control condition
    are used to train the noise filtration model, while the experimental
    communication scores are the target of the perturbation model. Large
    residuals between predicted and observed scores highlight
    significantly perturbed interactions.

2.  **Communication Pattern Detection via Tensor Decomposition:**  
    The pipeline embeds the filtered cell–cell communication residuals
    into a four-dimensional tensor (indexed by sending cell type,
    receiving cell type, ligand–receptor pair, and condition). Canonical
    polyadic (CP) decomposition is then applied (again using Python) to
    extract latent factors that represent distinct, disease-specific
    signaling patterns.

3.  **Knowledge Graph–Driven Drug Target Discovery:**  
    After tensor decomposition, Python is used to run the MultiXRank
    algorithm—a random walk with restart—to propagate factor signals
    through a multiplex biological knowledge graph. Only after this step
    are the ChEMBL and MONDO identifiers updated to their respective
    drug and disease names in R.

This vignette details the functions and workflow within TrokaChatML,
highlighting the interplay between R and Python components that together
yield a robust and scalable analysis pipeline.

# 1. Folder Structure Setup

The analysis begins by establishing a standardized directory structure.
The function `TrokaChat_folder_structure` creates a main folder
(“TrokaChat”) with designated subdirectories for CSV outputs, DEG
results, tensor data, and visualization outputs.

``` r
TrokaChat_folder_structure("./")
```

# 2. Differential Expression and Perturbed Communication Analysis

## 2.1 Initial DEG Analysis with TrokaChat.DEG

TrokaChat.DEG performs a unique differential expression (DE) analysis
that extends beyond simply running Seurat’s `FindMarkers`. It leverages
species-specific ligand–receptor databases (by loading and filtering an
Excel file of signaling genes) and implements a two-dimensional DE
test: - **Differential Expression Across Condition and Cell Type:**  
For each control cluster $c_j$ in the control (c) sample, gene
expression is compared against all other control clusters (excluding
$c_j$) to identify genes enriched in that cluster. Similarly, each
experimental/diseased cluster $E_j$ is compared against all control
clusters (excluding $c_j$) to capture condition-driven changes. The
resulting average log2 fold-change (log2FC) values indicate whether a
gene is up- or downregulated in a specific cluster relative to the
control. - **Iterative DEG Analysis:**  
An initial round of DEG tests is performed for all valid cluster
comparisons, followed by a second round that aggregates and refines the
differential expression results. This iterative approach enables
detection of subtle, condition-specific expression differences that
standard methods might overlook.

``` r
TrokaChat.DEG(object = immune_cells,
              samples = 2,
              shortidents = c("NG", "DIAB"),
              control_condition = "NG",
              filepath = "./TrokaChat/",
              export_folder_1 = "csv1",
              clusters = "seurat_clusters",
              sample_clus = "sample_clus",
              cluster_range = 0:20,
              sample_species = "mouse",
              cluster_pct_thresh = 0,
              DefaultAssay = "SCT")
```

## 2.2 Perturbation Analysis via Python

To rigorously filter out technical noise and identify truly perturbed
cell–cell communications, TrokaChatML uses two LightGBM models executed
in Python: - **Noise Filtration Model:**  
A permutation-based approach is used to generate a null distribution of
communication scores by shuffling cell-type labels. The communication
scores computed for these permuted (null) datasets serve as training
data for the noise LightGBM model. - **ML-Based Disease Perturbation
Model:**  
The model is trained solely on control (NL) data (using a 90:10 split)
with the same set of categorical features (ligand, receptor, source,
target, pathway label, and composite features) to learn baseline
communication patterns. When applied to experimental (LS) data, the
model predicts expected communication scores; large deviations
(residuals) between the observed and predicted scores indicate
significantly perturbed interactions.

The function `perform_communication_analysis` then calculates the
percentage of cells expressing each gene within clusters and refines the
list of significant interactions based on these perturbation scores
computed in Python.

``` r
TC_overall_list <- readRDS(file = paste0("./TrokaChat/csv1/TC_overall_list.rds"))
perform_communication_analysis(seurat_obj = immune_cells, 
                               overall_list = TC_overall_list, 
                               n_permutations = 1000, 
                               output_dir = "./TrokaChat/", 
                               clusters_col = "seurat_clusters", 
                               sample_col = "sample")
```

## 2.3 Null Distribution Generation and Communication Score Computation

In R, the function `TrokaChat.DEG.nulldist` generates a null
distribution of communication scores by: - Creating a counts table of
cell numbers per cluster and condition. - Permuting cluster labels
multiple times (as specified by `n.perms`) to compute gene expression
percentages for each permutation. - Aggregating these permutation
results into CSV and HDF5 files that capture the expected distribution
of communication scores under random conditions.

The **Communication Score** is computed in Python using the following
formulation:

$$
\text{Communication Score} = \left( \text{log2FC}_L \times \text{pct.1}_L \times \frac{N_s}{N_{\text{total}}} \right) \times \left( \frac{1}{n} \sum_{i=1}^{n} \left( \text{log2FC}_{S_i} \times \text{pct.1}_{S_i} \times \frac{N_t}{N_{\text{total}}} \right) \right)
$$

where: - **log2FC**$_L$ and **pct.1**$_L$ are the ligand’s average log2
fold-change and the fraction of source cells expressing the ligand. -
**log2FC**$_{S_i}$ and **pct.1**$_{S_i}$ are the analogous values for
receptor subunit $S_i$ in the target cells. - $N_s$, $N_t$, and
$N_{\text{total}}$ denote the number of cells in the source cluster,
target cluster, and total number of cells, respectively.

This score integrates both the magnitude and prevalence of expression
changes to quantify potential cell–cell communication. In Python,
communication scores are computed for both the null distribution data
and the actual experimental data. The null distribution scores are used
to train the noise LightGBM model. Large residuals—i.e., significant
differences between the observed and predicted communication
scores—indicate interactions that are unlikely to be due to noise. In
the Perturbation model, the experimental scores are the target variable.
Large residuals—i.e., significant differences between the observed and
predicted communication scores—indicate likely perturbed interactions.

``` r
TrokaChat.DEG.nulldist(object = immune_cells, 
                       shortidents = c("NG", "DIAB"),
                       control_condition = "NG",
                       filepath = "./TrokaChat/",
                       export_folder = "csv1_test",
                       clusters = "seurat_clusters",
                       n.perms = 2,
                       assay = "SCT")
```

# 3. Tensor Construction and Decomposition

Module 2 of the pipeline organizes the filtered communication residuals
into a four-dimensional tensor. The dimensions represent: - **Sending
Cell Type** - **Receiving Cell Type** - **Ligand–Receptor Pair** -
**Condition (Control vs. Disease)**

The function `process_and_create_tensor` constructs this tensor,
computes summary statistics (e.g., the percentage of zero entries), and
exports both the tensor data and its dimensions for downstream analysis.
Tensor decomposition (using CP decomposition) is performed externally in
Python.

``` r
mapping <- create_cluster_mapping_from_seurat(immune_cells, 
                                              cluster_num_col = "seurat_clusters", 
                                              cluster_name_col = "cluster_names")
import_dir <- "./TrokaChat/cell_cell comm final/exp_final_with_residuals.csv"
prepared_data <- prepare_data(import_dir, mapping)
results <- process_and_create_tensor(result = prepared_data,
                                     tensor_data_file = "./TrokaChat/TEST/2condition_tensor_data.csv",
                                     tensor_dim_file = "./TrokaChat/TEST/2condition_tensor_dimensions.csv")
tensor <- results$tensor
mode_names <- results$mode_names
```

# 4. Factor Mapping and Visualization

After tensor decomposition is complete in Python, the next step is to
import the factor matrices into R. The function `load_and_map_factors`
loads these matrices and maps their indices to biologically meaningful
categories. Then, `generate_and_export_heatmaps` creates
publication-quality heatmaps that display the latent factors—with custom
labels and dimensions—offering visual insight into the disease-specific
communication patterns.

``` r
result <- load_and_map_factors(dir_path = "./TrokaChat/Decomposed Tensor",
                               num_modes = 4,
                               factor_prefix = "factor_",
                               mapping_prefix = "mapping_")
factors_mapped <- result$factors_mapped

# Define mode labels and custom heatmap sizes.
mode_labels <- c(
  "Source" = "Sending Cluster",
  "Target" = "Receiving Cluster",
  "LigandReceptor" = "Ligand-Receptor Pair",
  "Tools" = "Analytical Tool"
)
custom_sizes <- list(
  Source = c(10, 8),
  Target = c(10, 8),
  LigandReceptor = c(12, 10),
  Tools = c(8, 6)
)

generate_and_export_heatmaps(factors_mapped = factors_mapped,
                             mode_names = mode_names,
                             mode_labels = mode_labels,
                             cluster_mapping = mapping,
                             output_dir = "./TrokaChat/TEST/Heatmaps/",
                             heatmap_sizes = custom_sizes)
```

# 5. Knowledge Graph–Driven Drug Target Discovery

Once tensor decomposition is complete, Python is used to run the
MultiXRank algorithm—a random walk with restart—to propagate factor
signals through a multiplex biological knowledge graph. This step
prioritizes candidate genes, pathways, and drug targets. After
MultiXRank has been executed in Python, the ChEMBL and MONDO identifiers
are updated in R via the functions `update_chembl_to_drug_names` and
`update_mondo_to_disease_names`.

Before these corrections, ligand–receptor factor processing and network
generation are performed. The function `process_ligand_receptor_factors`
refines the factor data, and `generate_factor_networks` converts these
data into network representations (monoplex and bipartite networks).

``` r
process_ligand_receptor_factors(factors_mapped = factors_mapped,
                                mode_names = mode_names,
                                species = "mouse",
                                output_dir = "./TrokaChat/TEST/Decomposed Tensor",
                                significant_data_file = "./TrokaChat/Decomposed Tensor/Cleaned Factors/cleaned_LigandReceptor_factors/significant_data.xlsx")

generate_factor_networks(excel_file_path = "./TrokaChat/Decomposed Tensor/Cleaned Factors/cleaned_LigandReceptor_factors/significant_data.xlsx",
                           species = "mouse",
                           kg_dir = "./TrokaChat/TEST/KG")
```

After MultiXRank has been run in Python, the external updates are
applied:

``` r
update_chembl_to_drug_names("/path/to/your/KG/TrokaChatML KG/", cores = 11)
update_mondo_to_disease_names("/path/to/your/KG/TrokaChat/TEST/KG/TrokaChatML KG/")
```

# 6. Integration with Python

A key strength of TrokaChatML is its hybrid workflow. The LightGBM
models for noise filtration and perturbation analysis are executed in
Python before the filtered communications are imported into R for tensor
assembly. Similarly, after the tensor decomposition, Python runs the
MultiXRank algorithm to perform the random walk with restart on the
multiplex knowledge graph. Only after these Python-based steps are
complete does the pipeline return to R for visualization and final
identifier updates (e.g., mapping ChEMBL and MONDO IDs to drug and
disease names).

*Note:* Detailed instructions for running the Python components are
provided in the repository’s README.

# 7. Running the Complete Pipeline

Below is a consolidated script that illustrates how to run the entire
TrokaChatML workflow from data loading to final output generation.

``` r
# Load necessary libraries and your scRNAseq data.
immune_cells <- readRDS("./immune_cells.rds")

# Visualize the initial data.
DimPlot(immune_cells, split.by = "sample")
Idents(immune_cells) <- "seurat_clusters"
FeaturePlot(immune_cells, split.by = "sample", features = c("Tgfb1"), label = TRUE)

# Set up the folder structure.
TrokaChat_folder_structure("./")

# Preprocess data: subset samples and reformat metadata.
immune_cells <- subset(immune_cells, subset = sample %in% c("NG Cre minus", "DIAB Cre minus"))
immune_cells@meta.data <- immune_cells@meta.data %>% mutate(sample = gsub(" Cre minus", "", sample))
immune_cells$sample_clus <- paste(immune_cells$sample, immune_cells$seurat_clusters, sep = "_")

# Run the initial DEG analysis.
TrokaChat.DEG(object = immune_cells,
              samples = 2,
              shortidents = c("NG", "DIAB"),
              control_condition = "NG",
              filepath = "./TrokaChat/",
              export_folder_1 = "csv1",
              clusters = "seurat_clusters",
              sample_clus = "sample_clus",
              cluster_range = 0:20,
              sample_species = "mouse",
              cluster_pct_thresh = 0,
              directory_path = "./TrokaChat/",
              DefaultAssay = "SCT")

# Perform communication analysis (with Python-based LightGBM models executed externally).
TC_overall_list <- readRDS(file = paste0("./TrokaChat/csv1/TC_overall_list.rds"))
perform_communication_analysis(seurat_obj = immune_cells, 
                               overall_list = TC_overall_list, 
                               n_permutations = 1000, 
                               output_dir = "./TrokaChat/", 
                               clusters_col = "seurat_clusters", 
                               sample_col = "sample")

# Generate the null distribution in R and transfer it to Python.
TrokaChat.DEG.nulldist(object = immune_cells, 
                       shortidents = c("NG", "DIAB"),
                       control_condition = "NG",
                       filepath = "./TrokaChat/",
                       export_folder = "csv1_test",
                       clusters = "seurat_clusters",
                       n.perms = 2,
                       assay = "SCT")

# Tensor construction.
mapping <- create_cluster_mapping_from_seurat(immune_cells, 
                                              cluster_num_col = "seurat_clusters", 
                                              cluster_name_col = "cluster_names")
import_dir <- "./TrokaChat/cell_cell comm final/exp_final_with_residuals.csv"
prepared_data <- prepare_data(import_dir, mapping)
results <- process_and_create_tensor(result = prepared_data,
                                     tensor_data_file = "./TrokaChat/TEST/2condition_tensor_data.csv",
                                     tensor_dim_file = "./TrokaChat/TEST/2condition_tensor_dimensions.csv")
tensor <- results$tensor
mode_names <- results$mode_names

# Factor mapping and heatmap generation.
result <- load_and_map_factors(dir_path = "./TrokaChat/Decomposed Tensor",
                               num_modes = 4,
                               factor_prefix = "factor_",
                               mapping_prefix = "mapping_")
factors_mapped <- result$factors_mapped

generate_and_export_heatmaps(factors_mapped = factors_mapped,
                             mode_names = mode_names,
                             mode_labels = mode_labels,
                             cluster_mapping = mapping,
                             output_dir = "./TrokaChat/TEST/Heatmaps/",
                             heatmap_sizes = custom_sizes)

# Ligand–receptor factor processing.
process_ligand_receptor_factors(factors_mapped = factors_mapped,
                                mode_names = mode_names,
                                species = "mouse",
                                output_dir = "./TrokaChat/TEST/Decomposed Tensor",
                                significant_data_file = "./TrokaChat/Decomposed Tensor/Cleaned Factors/cleaned_LigandReceptor_factors/significant_data.xlsx")

# Network generation.
generate_factor_networks(excel_file_path = "./TrokaChat/Decomposed Tensor/Cleaned Factors/cleaned_LigandReceptor_factors/significant_data.xlsx",
                           species = "mouse",
                           kg_dir = "./TrokaChat/TEST/KG")

# At this point, Python runs MultiXRank to execute the random walk with restart.
# Once complete, update external databases.
update_chembl_to_drug_names("/path/to/your/KG/TrokaChatML KG/", cores = 11)
update_mondo_to_disease_names("/path/to/your/KG/TrokaChat/TEST/KG/TrokaChatML KG/")
```

# Conclusion

TrokaChatML combines rigorous noise modeling, machine learning–based
perturbation analysis (with key steps executed in Python), tensor
decomposition, and multiplex knowledge graph integration into a unified
pipeline. By doing so, it not only recapitulates known cell–cell
communication patterns in diseases like atopic dermatitis but also
uncovers novel therapeutic targets—such as non-canonical NF-κB dynamics
in dendritic cells—that may have been overlooked by existing methods.
This modular, scalable approach provides a robust framework for both
mechanistic insight and drug target discovery.

------------------------------------------------------------------------

*End of Vignette*
