# Caire et al., 2026 (in press)

This directory contains metadata and scripts used in `Caire et al. 2026` (paper in press), required to replicate the analysis described in the paper.
<br/><br/>


## Contents

### Example script for the analysis of scRNA-seq samples (.R)


### Integration of scRNA-seq samples (.R)


### Differential expression analysis of human tumor pseudobulks (.R)


### 2D morphology in TNBC H&E-derived tumor areas (.R)


### Metadata of scRNA-seq data (.csv) 


<br/><br/>

## Code reproducibility
2D morphology analysis in TNBC H&E-derived tumor areas was performed using R (v...) and sf package (v.1.0-21). Tumor areas were identified using a Stardist-based weakly supervised machine learning object classifier in QuPath (v.0.5.0).  

All the other scripts were performed in R (v.3.6.3), using Seurat (v.3.1.5), popsicleR (v.0.2.1), scrublet (v.0.2.1), SingleR (v.1.0.6), infercnv (v.1.2.1).  

For pseudobulk analysis Matrix.utils R package (v.0.9.8) was used to aggregate scRNA-seq data in pseudobulk, and edgeR (v.3.28.1) was used to compute differential expression analysis  

For more details refer to methods section of the paper. 

<br/><br/>


## Citation

`Caire et al., 2026` paper (in press):
**A 3D Morphogenetic Program for Metastatic Outgrowth in Breast Cancer**. 

Robin Caire, Roberta Bordo, Francesca Zanconato, Tito Panciera, Paolo Contessotto, Mikaela Fakiola, Ramona Bason, Oriana Romano, Ambela Suli, Giusy Battilana, Matteo Marchionni, Mattia Forcato, Estelle Audoux, Sara Donzelli, Maria Vittoria Dieci, Gaia Griguolo, Maria AntoniaCarosi, Matteo Fassan, Vincenza Guzzardo, Angelo Paolo Dei Tos, Silvia Marsoni, Pei-Hsun Wu, Denis Wirtz, Shanshan He, Cecilia Casali, Francesco Volpin, Giovanni Blandino, Claudio Tripodo, Silvio Bicciato, Valentina Guarneri, Massimiliano Pagani, Michelangelo Cordenonsi and Stefano
Piccolo 
