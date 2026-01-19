# SpaceTrooper Data

<img src="https://github.com/bicciatolab/bicciatolab_data/blob/main/SpaceTrooper_data_and_scripts/SpaceTrooper_logo.png" alt="SpaceTrooper Logo" width="100" style="float:right; margin:0 0 10px 10px;"/>

This directory contains metadata files, InSituType reference matrices, and R scripts required to replicate the analysis described in the `SpaceTrooper` bioRxiv paper:
**[https://doi.org/10.64898/2025.12.24.696336](https://doi.org/10.64898/2025.12.24.696336)**


These files support the developmental version of `SpaceTrooper` package, available on [GitHub](https://github.com/drighelli/SpaceTrooper).

<br/><br/>

## Contents

### Metadata Files (.rds)

The [Spe_metadata](https://github.com/bicciatolab/bicciatolab_data/tree/main/SpaceTrooper_data_and_scripts/Spe_metadata) subdirectory contains 5 R Data Serialization (.rds) files with metadata for different spatial transcriptomics datasets:

1. **CosMx_protein_tonsil_metadata.rds** (~245 MB) - Metadata for CosMx protein tonsil dataset
2. **CosMx_rna_DBKERO_metadata.rds** (~6 MB) - Metadata for CosMx RNA DBKERO dataset
3. **CosMx_rna_pancreas_metadata.rds** (~6 MB) - Metadata for CosMx RNA pancreas dataset
4. **MERFISH_mouse_liver_metadata.rds** (~75 MB) - Metadata for MERFISH mouse liver dataset
5. **Xenium_lung_cancer_metadata.rds** (~12 MB) - Metadata for Xenium lung cancer dataset

### InSituType Reference Matrices

The [InSituType_references](https://github.com/bicciatolab/bicciatolab_data/tree/main/SpaceTrooper_data_and_scripts/InSituType_references) subdirectory contains reference matrices for performing InSituType cell phenotyping:

1. **Human DCIS scRNAseq simple_profileMatrix.csv** (~191 KB) - Reference matrix for human datasets (matches with CosMx RNA DBKERO dataset)
2. **Mouse liver cells_profileMatrix.csv** (~7.7 MB) - Reference matrix for mouse liver dataset (matches with MERFISH mouse liver dataset)

### Scripts

The [Scripts](https://github.com/bicciatolab/bicciatolab_data/tree/main/SpaceTrooper_data_and_scripts/Scripts) subdirectory contains the R scripts used to build the **SpatialExperiment** (Spe) objects and create paper's figures:

1. **Scripts_for_spe** - subdirectory containing the scripts used to build the **SpatialExperiment** objects used in the paper's analysis
2. **Scripts_for_figures** - subdirectory containing the scripts used to create the paper's figures

To reproduce the analysis, we provide a [requirements.txt](https://github.com/bicciatolab/bicciatolab_data/tree/main/SpaceTrooper_data_and_scripts/Scripts/requirements.txt) file containing the version of used packages, along with a [SpaceTrooper.yml](https://github.com/bicciatolab/bicciatolab_data/tree/main/SpaceTrooper_data_and_scripts/Scripts/SpaceTrooper.yml) to create a reproducible `anaconda` environment using either **conda** or **mamba** bash command:

```bash
# conda env create -f SpaceTrooper.yml
mamba env create -f SpaceTrooper.yml
```

Downstream analysis of **figure 2** (and related supplementary images) was performed using Seurat package. A separate [requirements_Seurat.txt](https://github.com/bicciatolab/bicciatolab_data/tree/main/SpaceTrooper_data_and_scripts/Scripts/requirements_Seurat.txt) file is provided to reproduce this analysis, along with a separated  [Seurat_v5.yml](https://github.com/bicciatolab/bicciatolab_data/tree/main/SpaceTrooper_data_and_scripts/Scripts/Seurat_v5.yml):

```bash
# conda env create -f Seurat_v5.yml
mamba env create -f Seurat_v5.yml
```

<br/><br/>

## How to Download 

### For Files < 50 MB

Small files can be downloaded directly from GitHub:

```bash
# Clone the repository normally
git clone https://github.com/bicciatolab/bicciatolab_data.git
```

### For Large Files (> 50 MB)

Two files in this directory are stored using Git Large File Storage (LFS):
- `CosMx_protein_tonsil_metadata.rds` (~245 MB)
- `MERFISH_mouse_liver_metadata.rds` (~75 MB)

To download these large files, you need Git LFS installed:

#### 1. Install Git LFS

**On Ubuntu/Debian:**
```bash
sudo apt-get install git-lfs
```

**On macOS:**
```bash
brew install git-lfs
```

**On Windows:**
Download and install from [git-lfs.github.com](https://git-lfs.github.com)

#### 2. Initialize Git LFS
```bash
git lfs install
```

#### 3. Clone the repository with LFS files
```bash
git clone https://github.com/bicciatolab/bicciatolab_data.git
cd bicciatolab_data/SpaceTrooper_data
```

Git LFS will automatically download the large files during the clone process.

#### Alternative: Download LFS files after cloning
If you already cloned the repository without Git LFS:
```bash
cd bicciatolab_data
git lfs pull
```

<br/><br/>

## Citation

If you use these datasets, please cite the `SpaceTrooper` paper: [https://doi.org/10.64898/2025.12.24.696336](https://doi.org/10.64898/2025.12.24.696336)
