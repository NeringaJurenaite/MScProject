## MScProject
# Cell-type Specific Differential Expression Between Human Cancer Types

This code has been implemented using R and has been used to write a Masters Thesis on using multiple linear regression model with an interaction term to infer cell-type specific differential expression between Breast Cancer subtypes Luminal A, Luminal B, HER2-enriched, Basal-like/Triple-like, and Normal like.

# Required Packages 
Please install the following packages through R (version 4.0) prompt:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ggplot2", "EDASeq", "TCGAbiolinks", "dplyr", "DT", "SummarizedExperiment", "limma", "edgeR", "jtools", "MASS", "splines", "interactions", "readxl", "preprocessCore", "sjPlot", "DESeq2", "sjlabelled", "sjmisc", "NOISeq", "tidyverse", "org.Hs.eg.db", "dplyr", "gage", "fgsea", "GSVAdata", "GSEABase", "ggsci", "gridExtra", "ensembldb", "EnsDb.Hsapiens.v86", "xtable")
```

# Files
Expression Count Data was obtained from the TCGA database in HTSeq Count format and restricted due to patient privacy therefore, it is not uploaded here. It is in the format of a matrix where Ensembl gene names are rows and sample names are columns.
Sample subtype data and tumour purity data has been obtained using TCGAbiolinks package, but the tumour purity data can also be accessed form [this link](https://www.nature.com/articles/ncomms9971#Sec14)
