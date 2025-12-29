# SpaceTrooper Data

This directory contains metadata files required to replicate the analysis described in the bioRxiv paper:
**doi: [https://doi.org/10.64898/2025.12.24.696336](https://doi.org/10.64898/2025.12.24.696336)**

These files support the [SpaceTrooper](https://github.com/drighelli/SpaceTrooper) package analysis ([Bioconductor page](https://bioconductor.org/packages/3.23/bioc/html/SpaceTrooper.html)).

## Contents

### Metadata Files (.rds)

This directory contains 5 R Data Serialization (.rds) files with metadata for different spatial transcriptomics datasets:

1. **CosMx_protein_tonsil_metadata.rds** (~245 MB) - Metadata for CosMx protein tonsil dataset
2. **CosMx_rna_DBKERO_metadata.rds** (~6 MB) - Metadata for CosMx RNA DBKERO dataset
3. **CosMx_rna_pancreas_metadata.rds** (~6 MB) - Metadata for CosMx RNA pancreas dataset
4. **MERFISH_mouse_liver_metadata.rds** (~75 MB) - Metadata for MERFISH mouse liver dataset
5. **Xenium_lung_cancer_metadata.rds** (~12 MB) - Metadata for Xenium lung cancer dataset

### InSituType Reference Matrices

The `InSituType_references/` subdirectory contains reference matrices for performing InSituType cell phenotyping:

1. **Human DCIS scRNAseq simple_profileMatrix.csv** (~191 KB) - Reference matrix for human datasets (matches with Xenium lung cancer dataset)
2. **Mouse liver cells_profileMatrix.csv** (~7.7 MB) - Reference matrix for mouse liver dataset (matches with MERFISH mouse liver dataset)

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

## Usage

Load the metadata files in R:
```R
# Load metadata
metadata <- readRDS("SpaceTrooper_data/CosMx_protein_tonsil_metadata.rds")

# Load reference matrix for InSituType
reference <- read.csv("SpaceTrooper_data/InSituType_references/Mouse liver cells_profileMatrix.csv")
```

## Citation

If you use these datasets, please cite the SpaceTrooper paper:
- bioRxiv doi: [https://doi.org/10.64898/2025.12.24.696336](https://doi.org/10.64898/2025.12.24.696336)
