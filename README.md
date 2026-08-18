
# AML-Gene-Analysis

R code for the Multilayer Feature Screening and Clustering Algorithm for High-Dimensional Genomic Data, with Application to Acute Myeloid Leukemia (AML).

## Data Source
The microarray gene expression dataset analyzed in this project is publicly available from the NCBI Gene Expression Omnibus (GEO) under accession number [GSE9476](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9476).

## Analytical Workflow
The repository provides reproducible R scripts for the following steps:
1. **Data Preprocessing & Quality Control** (`S2_File_12_5_1405.R`)
2. **Supervised Feature Screening & Clustering** (`FDR1_svm_rf_Lasso.R`): Non-parametric Mann-Whitney U test, modified FDR control, and hierarchical clustering of candidate genes.
3. **Regularization & Binary Classification** (`Boruta_svm_rf_Lasso.R`): Feature selection and regularized classification models (LASSO, Elastic Net, SVM, Random Forest).

## Citation
If you use this code or methodology in your research, please cite our corresponding publication.
