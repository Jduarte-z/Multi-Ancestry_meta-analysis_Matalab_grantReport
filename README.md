# Multi-Ancestry_meta-analysis_Matalab_grantReport

This repository contains the methods employed to conduct a Trans-ethnic meta-regression of genome-wide association studies for Parkinson’s Disease risk. As part of a series of in house 

We used eight different GWAS summary statistics from ancestrally diverse datasets. And used MR-MEGA (https://pubmed.ncbi.nlm.nih.gov/28911207/) to perform a meta-analysis with meta-regression. 

Briefly, MR-MEGA performs a fixed-effects meta-analysis across the input summary statistics to get the pooled effect sizes for variants with a minor allele frequency >1%. Then, using the allele frequency information of each variant, builds a matrix of mean pairwise allele frequency differences between the datasets, performs a principal component analysis of that matrix and then runs a meta-regression. Regressing each variant’s effect sizes on the principal components (PCs) of that allele frequency matrix. Allowing the possibility of recovering the significance of variants that have heterogeneous effects across cohorts (where the heterogeneity is correlated with ancestry). 

The main steps undertaken consisted on: 
1) Formatting of summary statistics
2) Run the first round of MR-MEGA to decide how many PCs to include
3) Run the second round of MR-MEGA
4) Format the output summary statistics
5) Interrogate the discovery of novel loci
6) Plot findings 


## 1) Format input summary statistics 

The following scripts were used to format the different summary statistics intended to be used. 
Considering that each file has its own way of representing the data, some specific modifications were applied to each file. 

### LARGE-PD
Starting with the summary statistics from the GWAS conducted in LARGE-PD phase 1 and phase 2 (with the latest genetic data available to date). 

Modifications made to LARGE-PD sumstats: 
The SNPIDs for these two GWAS had to be modified. Adding the “chr” prefix. 
Add the sample size of each GWAS. 

<details>
    <summary>editPhase1.py</summary>
            
```python=

import pandas as pd
import numpy as np
import gzip

# === File paths ===
INPUT_FILE = "saige_phase1.tsv"  # or .gz #the zipped file should not output any excluded rows because the one that has missing values is the plain .txt file 
OUTPUT_FILE = "lpd_phase1_saige_modified.tsv"
EXCLUDED_ROWS_FILE = "excluded_rows_lpd_phase1_saige.txt"
EXCLUSION_SUMMARY_FILE = "exclusion_summary_lpd_phase1_saige.txt"

# === Load file with gzip support ===
if INPUT_FILE.endswith(".gz"):
    with gzip.open(INPUT_FILE, 'rt') as f:
        df = pd.read_csv(f, sep="\t")
else:
    df = pd.read_csv(INPUT_FILE, sep="\t")

# # === Add sample size column ===
df["N"] = 1478


marker = (
        "chr"
        + df["CHR"].astype("string")
        + ":"
        + df["POS"].astype("Int64").astype("string")
        + ":"
        + df["Allele1"]
        + ":"
        + df["Allele2"]
    )
df.insert(0, "MARKERNAME", marker)

# === Rename columns ===
df = df.rename(columns={
    "Allele2": "EA",
    "Allele1": "NEA",
    "AF_Allele2": "EAF",
    "CI95_LO": "OR_95L",
    "CI95_HI": "OR_95U"
})

# === Reorder columns ===
final_cols = ["MARKERNAME", "CHR", "POS", "N", "EA", "NEA", "EAF", "OR", "OR_95L", "OR_95U"]
df = df[final_cols]

# === Handle missing values ===
is_missing = df.isnull()
excluded_rows = df[is_missing.any(axis=1)].copy()
included_rows = df.dropna()

# === Log excluded rows with reasons ===
with open(EXCLUDED_ROWS_FILE, "w") as excl_file:
    excl_file.write("MARKERNAME\tMissing_Columns\n")
    for idx, row in excluded_rows.iterrows():
        missing_cols = is_missing.loc[idx]
        missing_fields = ",".join(missing_cols[missing_cols].index)
        excl_file.write(f"{row['MARKERNAME']}\t{missing_fields}\n")

# === Summary of missingness by column ===
missing_summary = is_missing.sum()
missing_summary = missing_summary[missing_summary > 0]
with open(EXCLUSION_SUMMARY_FILE, "w") as summary_file:
    summary_file.write("Column\tMissing_Entries\n")
    for col, count in missing_summary.items():
        summary_file.write(f"{col}\t{count}\n")

# === Save final cleaned output ===
included_rows.to_csv(OUTPUT_FILE, sep="\t", index=False)

```
</details>

<details>
    <summary>editPhase2.py</summary>
            
```python=

import pandas as pd
import numpy as np
import gzip

# === File paths ===
INPUT_FILE = "saige_phase2.tsv"  # or .gz #the zipped file should not output any excluded rows because the one that has missing values is the plain .txt file 
OUTPUT_FILE = "lpd_phase2_batch1-8_saige_modified.tsv"
EXCLUDED_ROWS_FILE = "excluded_rows_lpd_phase2_batch1-8_saige.txt"
EXCLUSION_SUMMARY_FILE = "exclusion_summary_lpd_phase2_batch1-8_saige.txt"

# === Load file with gzip support ===
if INPUT_FILE.endswith(".gz"):
    with gzip.open(INPUT_FILE, 'rt') as f:
        df = pd.read_csv(f, sep="\t")
else:
    df = pd.read_csv(INPUT_FILE, sep="\t")

# # === Add sample size column ===
df["N"] = 5143

marker = (
        "chr"
        + df["CHR"].astype("string")
        + ":"
        + df["POS"].astype("Int64").astype("string")
        + ":"
        + df["Allele1"]
        + ":"
        + df["Allele2"]
    )
df.insert(0, "MARKERNAME", marker)

# === Rename columns ===
df = df.rename(columns={
    "Allele2": "EA",
    "Allele1": "NEA",
    "AF_Allele2": "EAF",
    "CI95_LO": "OR_95L",
    "CI95_HI": "OR_95U"
})

# === Reorder columns ===
final_cols = ["MARKERNAME", "CHR", "POS", "N", "EA", "NEA", "EAF", "OR", "OR_95L", "OR_95U"]
df = df[final_cols]

# === Handle missing values ===
is_missing = df.isnull()
excluded_rows = df[is_missing.any(axis=1)].copy()
included_rows = df.dropna()

# === Log excluded rows with reasons ===
with open(EXCLUDED_ROWS_FILE, "w") as excl_file:
    excl_file.write("MARKERNAME\tMissing_Columns\n")
    for idx, row in excluded_rows.iterrows():
        missing_cols = is_missing.loc[idx]
        missing_fields = ",".join(missing_cols[missing_cols].index)
        excl_file.write(f"{row['MARKERNAME']}\t{missing_fields}\n")

# === Summary of missingness by column ===
missing_summary = is_missing.sum()
missing_summary = missing_summary[missing_summary > 0]
with open(EXCLUSION_SUMMARY_FILE, "w") as summary_file:
    summary_file.write("Column\tMissing_Entries\n")
    for col, count in missing_summary.items():
        summary_file.write(f"{col}\t{count}\n")

# === Save final cleaned output ===
included_rows.to_csv(OUTPUT_FILE, sep="\t", index=False)

```
</details>


### LARGE-PD
