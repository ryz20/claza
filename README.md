# claza_2025

Repository for analysis of n=58 human PBMC bulk RNA-sequencing samples and n=58 paired kidney microarray samples from Doberer et al. (2021). This is not a supported package.

# Analysis pipeline
This repository contains the following files which illustrate the methodology in our paper and accompanies supplementary material and sequencing data available on GEO and referenced in our paper.
The following files include, for PBMC and kidney samples respectively:
explore.R (key data exploration including clustering samples and PCA)
de_analysis.R (differential expression analysis)
gsea.R (pathway analysis)
wgcna.R (weighted gene co-expression network analysis, and subsequent module analysis)
wgcna_gsea.R (gene set enrichment analysis of modules)
wgcna_pres.R (module preservation analysis of wgcna modules)

# Genesets
This folder contains genesets used from external sources or other analysis for use in the analysis. 
# Additional files
This folder contains additional files relevant to the analyses above.
