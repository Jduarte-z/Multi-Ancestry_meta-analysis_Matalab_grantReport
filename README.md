# Multi-Ancestry_meta-analysis_Matalab_grantReport

This repository contains the methods employed to conduct a Trans-ethnic meta-regression of genome-wide association studies for Parkinson’s Disease risk.

We used eight different GWAS summary statistics from ancestrally diverse datasets. And used MR-MEGA (https://pubmed.ncbi.nlm.nih.gov/28911207/) to perform a meta-analysis with meta-regression. 

Briefly, MR-MEGA performs a fixed-effects meta-analysis across the input summary statistics to get the pooled effect sizes for variants with a minor allele frequency >1%. Then, using the allele frequency information of each variant, builds a matrix of mean pairwise allele frequency differences between the datasets, performs a principal component analysis of that matrix and then runs a meta-regression. Regressing each variant’s effect sizes on the principal components (PCs) of that allele frequency matrix. Allowing the possibility of recovering the significance of variants that have heterogeneous effects across cohorts (where the heterogeneity is correlated with ancestry). 

The main steps undertaken consisted on: 
1) Formatting of summary statistics
2) Run the first round of MR-MEGA to decide how many PCs to include
3) Run the second round of MR-MEGA
4) Format the output summary statistics
5) Interrogate the discovery of novel loci
6) Plot findings 
