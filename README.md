# Multi-Ancestry_meta-analysis_Matalab_grantReport

This repository contains the methods employed to conduct a Trans-ethnic meta-regression of genome-wide association studies for Parkinson’s Disease risk. As part of a series of in house preliminar experiments.

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
Starting with the summary statistics from the GWAS conducted in LARGE-PD phase 1 and phase 2 (with the latest genetic data available to date). Both were in build hg38. 

For these summary statistics The SNPIDs had to be modified. Adding the “chr” prefix, and esuring a clean SNPID that consisted of chromosome, position, reference allele and alternative allele. 
Add the sample size of each GWAS. 

<details>
    <summary>editPhase1.py</summary>
            
```python=

import pandas as pd
import numpy as np
import gzip

# === File paths ===
INPUT_FILE = "saige_phase1.tsv"  # or .gz 
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
INPUT_FILE = "saige_phase2.tsv"  # or .gz 
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


### GP2 European Meta-analysis 

This dataset consists of the latest genome-wide association study in European subjects conducted by GP2. 

It was downloaded through this link: https://ndkp.hugeamp.org/research.html?pageid=a2f_downloads_280 

For this particular dataset, we had to compute the OR and 95%CI from the beta and standard error for each SNP, while adding the sample size of the study. 

<details>
    <summary>editGP2_eur.py</summary>
            
```python=

import pandas as pd
import numpy as np
import gzip

# === File paths ===
INPUT_FILE = "GP2_ALL_EUR_ALL_DATASET_HG38_12162024.txt"  # or .gz 
OUTPUT_FILE = "eurGP2_modified.tsv"
EXCLUDED_ROWS_FILE = "excluded_rows_eurGP2.txt"
EXCLUSION_SUMMARY_FILE = "exclusion_summary_eurGP2.txt"
EXCLUDED_ROWS_OVERFLOW = "exp_overflow_rows_eurGP2.txt"

# === Load file with gzip support ===
if INPUT_FILE.endswith(".gz"):
    with gzip.open(INPUT_FILE, 'rt') as f:
        df = pd.read_csv(f, sep="\t")
else:
    df = pd.read_csv(INPUT_FILE, sep="\t")

# # === Add sample size column ===
df["N"] = 1827641

b  = pd.to_numeric(df["beta(random)"], errors="coerce")
se = pd.to_numeric(df["standard_error"], errors="coerce")
z  = 1.96
hi = b + z*se
lo = b - z*se

def safe_exp(a):
    a = np.clip(a, -700, 700)  # stay below ~709.78 overflow; avoid tiny underflow too
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        return np.exp(a)

df["OR"]     = safe_exp(b)
df["OR_95U"] = safe_exp(hi)
df["OR_95L"] = safe_exp(lo)

overflow_mask = (b > 709) | (hi > 709)
if overflow_mask.any():
    df.loc[overflow_mask, ["SNP_ID","beta(random)","standard_error"]].to_csv(
        EXCLUDED_ROWS_OVERFLOW, sep="\t", index=False
    )



# === Rename columns ===
df = df.rename(columns={
    "SNP_ID": "MARKERNAME",
    "effect_allele": "EA",
    "other_allele": "NEA",
    "effect_allele_frequency": "EAF",
    "chromosome": "CHR",
    "base_pair_position": "POS"
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

### South African GWAS 

The South African summary statistics were shared directly from the authors of the pre-print: https://www.medrxiv.org/content/10.1101/2025.08.01.25331910v1.full.pdf+html 

For this particular dataset, we had to compute the OR and 95%CI. While updating each SNPID to match the format of chromosome, position, reference allele and alternative allele. 

<details>
    <summary>editSouthAfrica.py</summary>
            
```python=

import pandas as pd
import numpy as np
import gzip

# === File paths ===
INPUT_FILE = "Step_South_African_GWAS.txt"  # or .gz 
OUTPUT_FILE = "south_africanGWAS_modified.tsv"
EXCLUDED_ROWS_FILE = "excluded_rows_south_africanGWAS.txt"
EXCLUSION_SUMMARY_FILE = "exclusion_summary_south_africanGWAS.txt"
EXCLUDED_ROWS_OVERFLOW = "exp_overflow_rows_south_africanGWAS.txt"

# === Load file with gzip support ===
if INPUT_FILE.endswith(".gz"):
    with gzip.open(INPUT_FILE, 'rt') as f:
        df = pd.read_csv(f, sep="\t")
else:
    df = pd.read_csv(INPUT_FILE, sep="\t")

# # === Add sample size column ===
df["N"] = 1517


b  = pd.to_numeric(df["BETA"], errors="coerce")
se = pd.to_numeric(df["SE"], errors="coerce")
z  = 1.96
hi = b + z*se
lo = b - z*se

def safe_exp(a):
    a = np.clip(a, -700, 700)  # stay below ~709.78 overflow; avoid tiny underflow too
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        return np.exp(a)

df["OR"]     = safe_exp(b)
df["OR_95U"] = safe_exp(hi)
df["OR_95L"] = safe_exp(lo)

overflow_mask = (b > 709) | (hi > 709)
if overflow_mask.any():
    df.loc[overflow_mask, ["MarkerID","BETA","SE"]].to_csv(
        EXCLUDED_ROWS_OVERFLOW, sep="\t", index=False
    )

df["CHR"] = (
    df["CHR"].astype("string")
              .str.strip()
              .str.replace(r"^chr", "", case=False, regex=True)
)


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
    "AF_Allele2": "EAF"
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


### African GWAS 
