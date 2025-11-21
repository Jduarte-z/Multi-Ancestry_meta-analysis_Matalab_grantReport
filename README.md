# Multi-Ancestry_meta-analysis_Matalab_grantReport

This repository contains the methods employed to conduct a Trans-ethnic meta-regression of genome-wide association studies for Parkinson’s Disease risk. As part of a series of in house preliminar experiments.

The data repository accompaying this github could be found here:
https://1drv.ms/f/c/5E290E405B8F3AAB/IgClCBQoU0ntRK-S75UCccCDATqybV13aQOf4mkdflbzs-k?e=VVSDea

We used eight different GWAS summary statistics from ancestrally diverse datasets. And used MR-MEGA (https://pubmed.ncbi.nlm.nih.gov/28911207/) to perform a meta-analysis with meta-regression. 

Briefly, MR-MEGA performs a fixed-effects meta-analysis across the input summary statistics to get the pooled effect sizes for variants with a minor allele frequency >1%. Then, using the allele frequency information of each variant, builds a matrix of mean pairwise allele frequency differences between the datasets, performs a principal component analysis of that matrix and then runs a meta-regression. Regressing each variant’s effect sizes on the principal components (PCs) of that allele frequency matrix. Allowing the possibility of recovering the significance of variants that have heterogeneous effects across cohorts (where the heterogeneity is correlated with ancestry). 

The main steps undertaken consisted on: 
1) Formatting of summary statistics
2) Run the first round of MR-MEGA to decide how many PCs to include
3) Run the second round of MR-MEGA
4) Format the output summary statistics
5) Fine map genomic risk regions
6) Interrogate the discovery of novel loci
7) Plot findings 

Important disclaimers: 

This project was performed as part of the latest experiments that were supported by the NIH Grant R01 1R01NS112499-01A1.

The material outlined is preliminary and not intended for publication yet, since future multi-ancetry meta-analyses for PD risk will be conducted in conjunction the Global Parkinson’s Genetics Program (GP2).  

Regarding the analysis itself, further functionalities are under implementation, like alignment of the non-effect alleles to the reference allele in the reference genome, along with the performance of fixed and random effects meta-analysis. 

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
OUTPUT_FILE = "lpd_saige_phase2_modified.tsv"
EXCLUDED_ROWS_FILE = "excluded_rows_lpd_saige_phase2.txt"
EXCLUSION_SUMMARY_FILE = "exclusion_summary_lpd_saige_phase2.txt"

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
<details>
    <summary>fromat example</summary>
            
```python=

Format before:
MarkerID	CHR	POS	Allele1	Allele2	AF_Allele2	BETA	SE	OR	CI95_LO	CI95_HI	p.value
1:730869:C:T	1	730869	C	T	0.0198725	0.466893	0.483547	1.5950307259982397	0.6182520109112334	4.115025866440325	0.3342649
1:758443:G:C	1	758443	G	C	0.123096	-0.113392	0.174616	0.8928006136708357	0.6340424569430128	1.257160190208939	0.5160922

Format after:
MARKERNAME	CHR	POS	N	EA	NEA	EAF	OR	OR_95L	OR_95U
chr1:730869:C:T	1	730869	1478	T	C	0.0198725	1.5950307259982397	0.6182520109112334	4.115025866440325
chr1:758443:G:C	1	758443	1478	C	G	0.123096	0.8928006136708357	0.6340424569430128	1.257160190208939

```
</details>

### GP2 European Meta-analysis 

This dataset consists of the latest genome-wide association study in European subjects conducted by GP2. 

It was downloaded through this link: https://ndkp.hugeamp.org/research.html?pageid=a2f_downloads_280 

For this particular dataset, we had to compute the OR and 95%CI from the beta and standard error for each SNP, while adding the sample size of the study. 

They were in build hg38

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

<details>
    <summary>format example</summary>
            
```python=

Format before:
chromosome	base_pair_position	SNP_ID	effect_allele	other_allele	effect_allele_frequency	N_datasets	p_value	beta	standard_errop_value(random)	beta(random)	I
1	701203	chr1:701203:G:T	T	G	0.0107	2	0.2438	-0.1246	0.10690312944641893	0.7405	-0.0657	64.45
1	702040	chr1:702040:C:T	T	C	0.0504	2	0.5271	0.028	0.04427303845863247	0.5271	0.028	0.0


Format after:
MARKERNAME	CHR	POS	N	EA	NEA	EAF	OR	OR_95L	OR_95U
chr1:701203:G:T	1	701203	1827641	T	G	0.0107	0.9364117456810082	0.7593973401829615	1.1546879493126612
chr1:702040:C:T	1	702040	1827641	T	C	0.0504	1.028395684421425	0.9429187555345062	1.1216212187200558

```
</details>

### South African GWAS 

The South African summary statistics were shared directly from the authors of the pre-print: https://www.medrxiv.org/content/10.1101/2025.08.01.25331910v1.full.pdf+html 

For this particular dataset, we had to compute the OR and 95%CI. While updating each SNPID to match the format of chromosome, position, reference allele and alternative allele. 

They were in build hg 38

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

<details>
    <summary>format example</summary>
            
```python=

Format before:
CHR	POS	MarkerID	Allele1	Allele2	AC_Allele2	AF_Allele2	MissingRate	BETA	SE	Tstat	var	p.value	p.value.NA	Is.SPA	AF_case	AF_ctrl	N_case	N_ctrl	N_case_hom	N_case_het	N_ctrl_hom	N_ctrl_het
chr1	727233	rs151190501	G	A	45.637	0.0156184	0	0.00221495	0.510224	0.0085083	3.8413	9.965363E-01	9.965363E-01	false	0.0151561	0.0160198	679	782	0	18	0	21
chr1	730869	rs200188737	C	T	52.441	0.017947	0	-0.741933	0.457509	-3.54459	4.77751	1.048717E-01	1.048717E-01	false	0.0156922	0.0199047	679	782	0	16	1	21


Format after:
MARKERNAME	CHR	POS	N	EA	NEA	EAF	OR	OR_95L	OR_95U
chr1:727233:G:A	1	727233	1517	A	G	0.0156184	1.0022174048138466	0.36868078523639464	2.72441571878438
chr1:730869:C:T	1	730869	1517	T	C	0.017947	0.4761925451137218	0.19424196787696976	1.1674065213626248

```
</details>

### African GWAS 

This particular dataset was obtained through the NDKP portal: https://ndkp.hugeamp.org/dinspector.html?dataset=Rizig2023_Parkinsons_AF. 
It is important to note that this summary statistics do not contain the full original dtaset used in their publication. Since they not include 23andME
And we had to compute the OR and confidence intervals. 
They were in build hg38

<details>
    <summary>editAfricaNo23andMe.py</summary>
            
```python=

import pandas as pd
import numpy as np
import gzip

# === File paths ===
INPUT_FILE = "release5_11052023_summary_statistics_Rizig_et_al_2023_AFR_AAC_metaGWAS_no23andMe_hg38.txt"  # or .gz  
OUTPUT_FILE = "africaNo23me_modified.tsv"
EXCLUDED_ROWS_FILE = "excluded_rows_africaNo23me.txt"
EXCLUSION_SUMMARY_FILE = "exclusion_summary_africaNo23me.txt"
EXCLUDED_ROWS_OVERFLOW = "exp_overflow_rows_africaNo23me.txt"

# === Load file with gzip support ===
if INPUT_FILE.endswith(".gz"):
    with gzip.open(INPUT_FILE, 'rt') as f:
        df = pd.read_csv(f, sep="\t")
else:
    df = pd.read_csv(INPUT_FILE, sep="\t")

# # === Add sample size column ===
df["N"] = 3645


b  = pd.to_numeric(df["beta"], errors="coerce")
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
    df.loc[overflow_mask, ["variant_id","beta","standard_error"]].to_csv(
        EXCLUDED_ROWS_OVERFLOW, sep="\t", index=False
    )



# === Rename columns ===
df = df.rename(columns={
    "variant_id": "MARKERNAME",
    "effect_allele": "EA",
    "other_allele": "NEA",
    "effect_allele_frequency": "EAF",
    "chromosome": "CHR",
    "base_pair_location": "POS"
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
    <summary>format example</summary>
            
```python=

Format before
chromosome	base_pair_location	effect_allele	other_allele	beta	standard_error	effect_allele_frequency	p_value	variant_id	ref_allele	direction	HetISq	HetChiSq	HetDf	HetPVal	rsid
1	66861	T	C	-0.1072	0.4023	0.0724	0.7899	chr1:66861:C:T	C	??-	0.0	0.0	0	1.0	rs28375825
1	80346	C	G	0.4608	0.3456	0.8338	0.1824	chr1:80346:C:G	C	??+	0.0	0.0	0	1.0	rs376665626

Format after
MARKERNAME	CHR	POS	N	EA	NEA	EAF	OR	OR_95L	OR_95U
chr1:66861:C:T	1	66861	3645	T	C	0.0724	0.8983459858250229	0.4083184068469859	1.976461253558816
chr1:80346:C:G	1	80346	3645	C	G	0.8338	1.585341751221354	0.805271331975759	3.1210703378687366
```
</details>


### East Asia 

This dataset was obtained contacting the author of this paper directly: https://jamanetwork.com/journals/jamaneurology/fullarticle/2764340 
However, the initial summary statistics did not contain the allele frequencies per each SNP. Hence, very kindly, the author shared with us the allele frequencies of each of the east asian datasets that were meta-analyzed by them in order to compute the weighted allele frequencies. 

The original dataset was in build hg19

First we started editing the original sumstats computing the OR and confidence intervals, while making a brand new SNPID that follows the logic of chromosome, position, reference allele and alternative allele.  

<details>
    <summary>edit_foo.py</summary>
            
```python=

import pandas as pd
import numpy as np
import gzip

# === File paths ===
INPUT_FILE = "6724PDcases-24851controls-5843213snps-summary-stats-metaP-recalSE.txt"  # or .gz 
OUTPUT_FILE = "fooAll_modified.tsv"
EXCLUDED_ROWS_FILE = "excluded_rows_fooAll.txt"
EXCLUSION_SUMMARY_FILE = "exclusion_summary_fooAll.txt"
EXCLUDED_ROWS_OVERFLOW = "exp_overflow_rows_fooAll.txt"

# === Load file with gzip support ===
if INPUT_FILE.endswith(".gz"):
    with gzip.open(INPUT_FILE, 'rt') as f:
        df = pd.read_csv(f, sep="\t")
else:
    df = pd.read_csv(INPUT_FILE, sep="\t")

# # === Add sample size column ===
df["N"] = 31575

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
    df.loc[overflow_mask, ["SNP","BETA","SE"]].to_csv(
        EXCLUDED_ROWS_OVERFLOW, sep="\t", index=False
    )

marker = (
        "chr"
        + df["CHR"].astype("string")
        + ":"
        + df["BP"].astype("Int64").astype("string")
        + ":"
        + df["A2"]
        + ":"
        + df["A1"]
    )
df.insert(0, "MARKERNAME", marker)

# === Rename columns ===
df = df.rename(columns={
    "A1": "EA",
    "A2": "NEA",
    "BP": "POS"
})

# === Reorder columns ===
final_cols = ["MARKERNAME", "CHR", "POS", "N", "EA", "NEA", "OR", "OR_95L", "OR_95U"]
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
    <summary>format example</summary>
            
```python=
Format before:
CHR	BP	SNP	A1	A2	BETA	P	SE
1	794332	rs12127425	A	G	0.0013	0.9691	0.0335595656221841
1	832359	rs9697294:832359:C:T	T	C	-0.0932	0.1868	0.070600450487048

Format after
MARKERNAME	CHR	POS	N	EA	NEA	OR	OR_95L	OR_95U
chr1:794332:G:A	1	794332	31575	A	G	1.0013008453662857	0.9375579135244234	1.0693775482757104
chr1:832359:C:T	1	832359	31575	T	C	0.9110112798187364	0.7932817033173323	1.0462129007719911


```
</details>

After the initial formatting is done, then we had to append the weighted allele frequencies that were shared in different files corresponding to each east asian cohort meta-analyzed.  

<details>
    <summary>get_allele_frequencies.py</summary>
            
```python=

import pandas as pd

# --- Load base (the GWAS you want to append pooled EAF to) ---
fooAll = pd.read_csv("fooAll_modified.tsv", sep="\t", dtype={"MARKERNAME":"string"})
fooAll["MARKERNAME"] = fooAll["MARKERNAME"].str.strip()
for c in ("EA","NEA"):
    fooAll[c] = fooAll[c].str.upper().str.strip()

# --- Load cohorts (assumes each has MARKERNAME, EA, NEA, EAF, N) ---
fooChina   = pd.read_csv("../china/fooChina_modified.tsv", sep="\t", dtype={"MARKERNAME":"string"})
fooHk      = pd.read_csv("../hongKong/fooHK_modified.tsv", sep="\t", dtype={"MARKERNAME":"string"})
fooKorea   = pd.read_csv("../Korea/fooKorea_modified.tsv", sep="\t", dtype={"MARKERNAME":"string"})
fooSing    = pd.read_csv("../singapore/fooSingapore_modified.tsv", sep="\t", dtype={"MARKERNAME":"string"})
fooTaiwan  = pd.read_csv("../taiwan/fooTaiwan_modified.tsv", sep="\t", dtype={"MARKERNAME":"string"})

def prep(df, cohort_name):
    df = df.copy()
    df["MARKERNAME"] = df["MARKERNAME"].str.strip()
    for c in ("EA","NEA"):
        if c in df:
            df[c] = df[c].astype(str).str.upper().str.strip()
    # Keep only the needed columns; rename to generic names
    needed = {"MARKERNAME":"MARKERNAME", "EA":"EA_pop", "NEA":"NEA_pop", "EAF":"EAF_pop", "N":"N_pop"}
    missing = [k for k in needed if k not in df.columns]
    if missing:
        raise ValueError(f"{cohort_name}: missing columns {missing}. Rename or map them before running.")
    df = df[list(needed)].rename(columns=needed)
    df["cohort"] = cohort_name
    return df

cohorts = [
    prep(fooChina,  "china"),
    prep(fooHk,     "hong_kong"),
    prep(fooKorea,  "korea"),
    prep(fooSing,   "singapore"),
    prep(fooTaiwan, "taiwan"),
]

# --- Long table of all cohort contributions ---
long = pd.concat(cohorts, ignore_index=True)

# Merge with fooAll to access its EA/NEA for alignment
base = fooAll[["MARKERNAME", "EA", "NEA"]]
merged = base.merge(long, on="MARKERNAME", how="left")

# --- Harmonize EAF so it matches fooAll.EA ---
# If cohort's EA==fooAll.EA: keep EAF
# If cohort's EA==fooAll.NEA: flip (1-EAF)
# Try strand complements too (A<->T, C<->G), otherwise drop
COMPL = {"A":"T", "T":"A", "C":"G", "G":"C"}

def align_row(r):
    ea, nea = r["EA"], r["NEA"]
    pea, pnea, eaf = r["EA_pop"], r["NEA_pop"], r["EAF_pop"]
    if pd.isna(eaf) or pd.isna(pea) or pd.isna(pnea):
        return pd.NA

    # exact
    if pea == ea and pnea == nea:
        return eaf
    if pea == nea and pnea == ea:
        return 1 - eaf

    # complement
    cea  = COMPL.get(pea, pea)
    cnea = COMPL.get(pnea, pnea)
    if cea == ea and cnea == nea:
        return eaf
    if cea == nea and cnea == ea:
        return 1 - eaf

    return pd.NA

merged["EAF_aligned"] = merged.apply(align_row, axis=1)

# keep valid contributions
valid = merged[
    merged["EAF_aligned"].notna() & merged["N_pop"].notna() & (merged["N_pop"] > 0)
].copy()

# clamp to [0,1] defensively
valid["EAF_aligned"] = valid["EAF_aligned"].clip(0, 1)

# --- Weighted average per SNP ---
valid["wx"] = valid["N_pop"] * valid["EAF_aligned"]
agg = (
    valid.groupby("MARKERNAME")
         .agg(EAF=("wx", lambda s: s.sum() / valid.loc[s.index, "N_pop"].sum()),
              EAF_N_total=("N_pop", "sum"),
              EAF_sources=("cohort", lambda s: ",".join(sorted(s.unique()))))
         .reset_index()
)

# --- Append EAF to fooAll ---
fooAll = fooAll.merge(agg, on="MARKERNAME", how="left")

# Optional: basic QC counts
n_total   = fooAll.shape[0]
n_with    = fooAll["EAF"].notna().sum()
n_missing = n_total - n_with
print(f"Pooled EAF assigned to {n_with}/{n_total} SNPs; {n_missing} missing (no alignable cohort data).")

# Save
fooAll.to_csv("fooAll_with_EAF.tsv", sep="\t", index=False)

```
</details>

<details>
    <summary>format example with AF</summary>
            
```python=
Format before:
MARKERNAME	CHR	POS	N	EA	NEA	OR	OR_95L	OR_95U
chr1:794332:G:A	1	794332	31575	A	G	1.0013008453662857	0.9375579135244234	1.0693775482757104
chr1:832359:C:T	1	832359	31575	T	C	0.9110112798187364	0.7932817033173323	1.0462129007719911

Format after
MARKERNAME	CHR	POS	N	EA	NEA	OR	OR_95L	OR_95U	EAF	EAF_N_total	EAF_sources
chr1:794332:G:A	1	794332	31575	A	G	1.0013008453662855	0.9375579135244234	1.0693775482757104	0.14321473963578782	31575	china,hong_kong,korea,singapore,taiwan
chr1:832359:C:T	1	832359	31575	T	C	0.9110112798187364	0.7932817033173323	1.0462129007719911	0.026807011474267617	31575	china,hong_kong,korea,singapore,taiwan

```
</details>

After appending the allele frequencies for each SNP, the next step was to lift over the coordinates to build hg38, using an in house pipeline that leverages a custom table built with dbSNP reference files in order to perform a 1-1 conversion of coordinates. 

For the process of lifting over the summary statistics, we droped the MARKERNAME column, and build a new one once the coordinates were mapped to build hg38 

<details>
    <summary>format example liftover</summary>
            
```python=
Format before liftover (note that the MARKERNAME columns was dropped before using the liftover script)
CHR	POS	N	EA	NEA	OR	OR_95L	OR_95U	EAF	EAF_N_total	EAF_sources
1	794332	31575	A	G	1.0013008453662855	0.9375579135244234	1.0693775482757104	0.14321473963578782	31575	china,hong_kong,korea,singapore,taiwan
1	832359	31575	T	C	0.9110112798187364	0.7932817033173323	1.0462129007719911	0.026807011474267617	31575	china,hong_kong,korea,singapore,taiwan


Format after liftover
CHR	POS	N	EA	NEA	OR	OR_95L	OR_95U	EAF	EAF_N_total	EAF_sources
1	858952	31575	A	G	1.0013008453662855	0.9375579135244234	1.0693775482757104	0.1432147396357878	31575	china,hong_kong,korea,singapore,taiwan
1	896979	31575	T	C	0.9110112798187364	0.7932817033173323	1.0462129007719911	0.0268070114742676	31575	china,hong_kong,korea,singapore,taiwan
```
</details>

<details>
    <summary>liftover script</summary>
            
```python=

# consider that the custom reference files needed look like this:

RSID	A1	A2	CHR	POS_HG38	POS_HG37
rs1570391677	T	A	1	10001	10001
rs1570391692	A	C	1	10002	10002

# the script mainly uses the chromosome and position columns, as it is right now, we do not use the RSID as a key to lift over.

#liftover script
import argparse
from glob import glob
from pathlib import Path
import sys
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Lift over GWAS sumstats positions with tab-delimited reference maps.")
    p.add_argument("--sumstats", required=True, help="Path to the TAB-delimited sumstats file.")
    p.add_argument("--ref-glob", required=True, help='Glob for 22 TAB-delimited refs, e.g. "ref/chr*.tsv" or ".tsv.gz".')
    p.add_argument("--chr-col", default="CHR", help="Chromosome column name in sumstats (default: CHR).")
    p.add_argument("--pos-col", default="POS", help="Position column name in sumstats (default: POS).")
    p.add_argument("--from-build", choices=["hg37", "hg38"], required=True, help="Build of input coordinates.")
    p.add_argument("--to-build", choices=["hg37", "hg38"], required=True, help="Build to lift to.")
    p.add_argument("--out", required=True, help="Output TSV (mapped rows only). .gz supported via extension.")
    p.add_argument("--out-unmapped", required=True, help="Output TSV for unmapped rows. .gz supported via extension.")
    return p.parse_args()


def normalize_chr(s):
    # Remove optional chr prefix; coerce to pandas nullable Int64 so empties become <NA>.
    s = s.astype(str).str.replace(r"^chr", "", case=False, regex=True)
    s = pd.to_numeric(s, errors="coerce").astype("Int64")
    # Keep autosomes 1..22; others (e.g., X/Y/MT) become <NA> and won't map with autosome-only refs.
    return s.where(s.between(1, 22), other=pd.NA)


def to_int64_pos(s, name_for_errors):
    # Accept ints/floats/strings; blanks parse as NaN then become <NA>.
    if pd.api.types.is_integer_dtype(s):
        return s.astype("Int64")
    if pd.api.types.is_float_dtype(s):
        frac = s.dropna() % 1
        if not frac.eq(0).all():
            bad = s[~s.isna() & ~((s % 1) == 0)]
            raise ValueError(f"{name_for_errors}: found non-integer float positions, e.g. {bad.iloc[0]}")
        return s.round().astype("Int64")
    return pd.to_numeric(s, errors="coerce").astype("Int64")


def load_reference_map(ref_glob, from_build, to_build):
    # Read all TAB-delimited reference files and build a (CHR, POS_from) → POS_to series.
    files = sorted(glob(ref_glob))
    if not files:
        raise FileNotFoundError(f"No reference files matched: {ref_glob}")
    usecols = ["CHR", "POS_HG38", "POS_HG37"]
    frames = []
    for f in files:
        df = pd.read_csv(f, sep="\t", usecols=usecols, dtype=str, engine="c")
        df["CHR"] = normalize_chr(df["CHR"])
        df["POS_HG38"] = to_int64_pos(df["POS_HG38"], f"{f}::POS_HG38")
        df["POS_HG37"] = to_int64_pos(df["POS_HG37"], f"{f}::POS_HG37")
        frames.append(df)
    ref = pd.concat(frames, ignore_index=True)

    pos_from = "POS_HG38" if from_build == "hg38" else "POS_HG37"  # source column
    pos_to   = "POS_HG37" if from_build == "hg38" else "POS_HG38"  # target column

    # Keep only rows that can map (valid CHR + both source & target positions).
    ref = ref.dropna(subset=["CHR", pos_from, pos_to])

    # Curate duplicates: for each (CHR, POS_from), keep only if all target positions are identical.
    # Groups with conflicting targets are ambiguous and excluded.
    gcounts = ref.groupby(["CHR", pos_from])[pos_to].nunique(dropna=True)
    consistent_keys = gcounts[gcounts == 1].index  # unique mapping
    ambiguous_keys = gcounts[gcounts > 1].index    # conflicting mapping

    # Log a concise curation summary.
    print(
        f"[liftover] reference curation: total source keys={len(gcounts)}, "
        f"kept unique={len(consistent_keys)}, dropped ambiguous={len(ambiguous_keys)}",
        file=sys.stderr,
    )

    # Build a mapping Series; drop duplicate rows for identical keys just in case.
    ref_curated = ref.set_index(["CHR", pos_from]).loc[consistent_keys]
    ref_curated = ref_curated[~ref_curated.index.duplicated(keep="first")]
    ref_map = ref_curated[pos_to]

    # Sanity
    if ref_map.index.has_duplicates:
        raise ValueError("Unexpected duplicate keys remained after curation; please inspect reference files.")

    return ref_map.rename("_POS_TO_REF"), pos_to


def main():
    args = parse_args()
    if args.from_build == args.to_build:
        raise SystemExit("--from-build and --to-build must differ.")

    # Read sumstats and capture original column order (for later reordering).
    ss = pd.read_csv(args.sumstats, sep="\t", engine="c")
    for col in (args.chr_col, args.pos_col):
        if col not in ss.columns:
            raise SystemExit(f"Column '{col}' not found in sumstats.")
    original_cols = list(ss.columns)  # for column reordering in mapped

    # Build normalized join keys without altering the original columns.
    ss = ss.copy()
    ss["_CHR_KEY"] = normalize_chr(ss[args.chr_col])
    ss["_POS_FROM_KEY"] = to_int64_pos(ss[args.pos_col], f"{args.sumstats}::{args.pos_col}")

    # Load curated reference map and merge (left).
    ref_map, pos_to_name = load_reference_map(args.ref_glob, args.from_build, args.to_build)
    merged = ss.merge(
        ref_map,
        left_on=["_CHR_KEY", "_POS_FROM_KEY"],
        right_index=True,
        how="left",
        copy=False,
        validate="m:1",
    )

    # Create the target-build position column (Int64). 
    target_suffix = "_HG37" if args.to_build == "hg37" else "_HG38"
    new_pos_col = f"{args.pos_col}{target_suffix}"
    merged[new_pos_col] = merged["_POS_TO_REF"].astype("Int64")

    # Split into mapped vs unmapped; drop helper columns in both.
    mapped_mask = merged[new_pos_col].notna()
    mapped = merged.loc[mapped_mask].drop(columns=["_CHR_KEY", "_POS_FROM_KEY", "_POS_TO_REF"])
    unmapped = merged.loc[~mapped_mask].drop(columns=["_CHR_KEY", "_POS_FROM_KEY", "_POS_TO_REF"])

    # Drop the original source position column in the mapped output.
    if args.pos_col in mapped.columns:
        mapped = mapped.drop(columns=[args.pos_col])

    # Reinsert the new_pos_col where the original POS column used to be.
    #    - Remove the old POS
    #      to avoid duplicates, then insert new_pos_col at the original index.
    try:
        pos_idx = original_cols.index(args.pos_col)
    except ValueError:
        pos_idx = None  # Shouldn't happen due to earlier check, but be safe.

    if pos_idx is not None:
        base_order = [c for c in original_cols if c not in (args.pos_col, new_pos_col)]
        desired = base_order.copy()
        # Insert the liftover column exactly where POS used to be.
        desired.insert(pos_idx, new_pos_col)
        # Include any extra columns not in the original sumstats (rare, but safe).
        extras = [c for c in mapped.columns if c not in desired]
        mapped = mapped[desired + extras]

    # Report mapping stats.
    total = len(merged)
    n_mapped = int(mapped_mask.sum())
    n_unmapped = total - n_mapped
    print(f"[liftover] mapped {n_mapped}/{total} ({n_mapped/total:.2%}); unmapped {n_unmapped}", file=sys.stderr)

    # Write outputs (TAB; .gz inferred). Mapped has only the target POS; unmapped keeps original POS.
    mapped.to_csv(Path(args.out), sep="\t", index=False, compression="infer")
    unmapped.to_csv(Path(args.out_unmapped), sep="\t", index=False, compression="infer")
    print(f"[liftover] wrote mapped → {args.out}", file=sys.stderr)
    print(f"[liftover] wrote unmapped → {args.out_unmapped}", file=sys.stderr)


if __name__ == "__main__":
    main()

```
</details>

After the liftover process was done, we had to build from scratch the MARKERNAME to match with the appropiate coordinates. 

<details>
    <summary>edit_foo_afterLiftover</summary>
            
```python=
import pandas as pd
import numpy as np
import gzip

# === File paths ===
INPUT_FILE = "fooAll_finalCols_38.tsv"  # or .gz 
OUTPUT_FILE = "fooAll_hg38_modified_final.tsv"
EXCLUDED_ROWS_FILE = "excluded_rows_fooAll_hg38_modified_final.txt"
EXCLUSION_SUMMARY_FILE = "exclusion_summary_fooAll_hg38_modified_final.txt"

# === Load file with gzip support ===
if INPUT_FILE.endswith(".gz"):
    with gzip.open(INPUT_FILE, 'rt') as f:
        df = pd.read_csv(f, sep="\t")
else:
    df = pd.read_csv(INPUT_FILE, sep="\t")

# # === Add sample size column ===

marker = (
        "chr"
        + df["CHR"].astype("string")
        + ":"
        + df["POS"].astype("Int64").astype("string")
        + ":"
        + df["NEA"]
        + ":"
        + df["EA"]
    )
df.insert(0, "MARKERNAME", marker)



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
    <summary>format example</summary>
            
```python=
Format before
CHR	POS	N	EA	NEA	OR	OR_95L	OR_95U	EAF	EAF_N_total	EAF_sources
1	858952	31575	A	G	1.0013008453662855	0.9375579135244234	1.0693775482757104	0.1432147396357878	31575	china,hong_kong,korea,singapore,taiwan
1	896979	31575	T	C	0.9110112798187364	0.7932817033173323	1.0462129007719911	0.0268070114742676	31575	china,hong_kong,korea,singapore,taiwan

Format after
MARKERNAME	CHR	POS	N	EA	NEA	EAF	OR	OR_95L	OR_95U
chr1:858952:G:A	1	858952	31575	A	G	0.1432147396357878	1.0013008453662855	0.9375579135244234	1.0693775482757104
chr1:896979:C:T	1	896979	31575	T	C	0.0268070114742676	0.9110112798187364	0.7932817033173323	1.0462129007719911


```
</details>


### South Asia 

For South Asian populations we included two studies. 

The first one could be downloaded through the NDKP portal (Kishore et al., 2025): https://ndkp.hugeamp.org/research.html?pageid=a2f_downloads_280
And the second one can be downloaded through zenodo (Andrew et al., 2024): https://zenodo.org/records/8436983

Both datasets were in build hg38

For Kishore et al., dataset we had to compute the OR and confidence interval, while building a brand new SNPID 

<details>
    <summary>editKishore.py</summary>
            
```python=

import pandas as pd
import numpy as np
import gzip

# === File paths ===
INPUT_FILE = "PD_GWAS_India_4806Cases_6364Controls_Summary_Statistics.tsv"  # or .gz 
OUTPUT_FILE = "kishore_modified.tsv"
EXCLUDED_ROWS_FILE = "excluded_rows_kishore.txt"
EXCLUSION_SUMMARY_FILE = "exclusion_summary_kishore.txt"
EXCLUDED_ROWS_OVERFLOW = "exp_overflow_rows_kishore.txt"

# === Load file with gzip support ===
if INPUT_FILE.endswith(".gz"):
    with gzip.open(INPUT_FILE, 'rt') as f:
        df = pd.read_csv(f, sep="\t")
else:
    df = pd.read_csv(INPUT_FILE, sep="\t")

# # === Add sample size column ===
df["N"] = 11170

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
    df.loc[overflow_mask, ["SNP","BETA","SE"]].to_csv(
        EXCLUDED_ROWS_OVERFLOW, sep="\t", index=False
    )

marker = (
        "chr"
        + df["CHR"].astype("string")
        + ":"
        + df["POS"].astype("Int64").astype("string")
        + ":"
        + df["A2"]
        + ":"
        + df["A1"]
    )
df.insert(0, "MARKERNAME", marker)

# === Rename columns ===
df = df.rename(columns={
    "A1": "EA",
    "A2": "NEA",
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
    <summary>format example</summary>
            
```python=
Format before 
CHR	POS	SNP	A1	A2	EAF	BETA	SE	P
1	727242	rs61769339	A	G	0.181916	0.0326937	0.0353983	0.355696
1	727717	rs61769340	C	G	0.764772	-0.0208611	0.0322121	0.517232

Format after
chr1:727242:G:A	1	727242	11170	A	G	0.181916	1.0332340111898817	0.9639778440368273	1.107465828684282
chr1:727717:G:C	1	727717	11170	C	G	0.764772	0.9793549875302142	0.9194341835689046	1.0431809135889345

```
</details>

For Andrew et al., dataset, we made the following changes:

<details>
    <summary>editAndrew.py</summary>
            
```python=
import pandas as pd
import numpy as np
import gzip

# === File paths ===
INPUT_FILE = "PDDiagnosis_GWAS.txt"  # or .gz 
OUTPUT_FILE = "andrew_modified.tsv"
EXCLUDED_ROWS_FILE = "excluded_rows_adrew.txt"
EXCLUSION_SUMMARY_FILE = "exclusion_summary_adrew.txt"

# === Load file with gzip support ===
if INPUT_FILE.endswith(".gz"):
    with gzip.open(INPUT_FILE, 'rt') as f:
        df = pd.read_csv(f, sep="\t")
else:
    df = pd.read_csv(INPUT_FILE, sep=" ")

# # === Add sample size column ===
df["N"] = 1878

marker = (
        "chr"
        + df["CHR"].astype("string")
        + ":"
        + df["BP"].astype("Int64").astype("string")
        + ":"
        + df["A2"]
        + ":"
        + df["A1"]
    )
df.insert(0, "MARKERNAME", marker)

# === Rename columns ===
df = df.rename(columns={
    "A1": "EA",
    "A2": "NEA",
    "MAF": "EAF",
    "L95": "OR_95L",
    "U95": "OR_95U",
    "BP": "POS",
    "BETA": "OR" #note that the authors in their repository state that the BETA column contains the OR value with respect to A1
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
    <summary>Format example</summary>
            
```python=
Format before
CHR BP A1 A2 BETA SE L95 U95 STAT P MAF MAF.cases MAF.cont DR2
1 817341 A G 1.1822 0.0933979 0.984444 1.41968 1.7921 0.0731177 0.239354259850905 0.2612 0.2311 1
1 819123 G A 1.19381 0.0928384 0.9952 1.43205 1.90813 0.056374 0.247344568690096 0.2689 0.2392 1

Format after
MARKERNAME	CHR	POS	N	EA	NEA	EAF	OR	OR_95L	OR_95U
chr1:817341:G:A	1	817341	1878	A	G	0.239354259850905	1.1822	0.984444	1.41968
chr1:819123:A:G	1	819123	1878	G	A	0.247344568690096	1.19381	0.9952	1.43205


```
</details>


## 2) MR-MEGA first round run 

This step is meant to analyze the principal components (PCs) derived from the pairwise allele frequency matrix (a proxy for ancestry) that MR-MEGA computes for the meta-regression. And choose the optimal number of PCs to be included to maximize discovery power and cohort separation–too few PCs could give a hard time separating different datasets in terms of ancestry and too many could overfit the meta-regression compromising statistical power–. 

By default, MR-MEGA restricts the number of principal components to be computed in the following way: T < K-2. Where T is the number of PCs and K the number of studies. So, when working with eight studies, the default number of PCs is five. And these are the ones we plot in order to examine if they are the right amount. 

This step borrows the plotting scripts from the latest multi-macentry meta-analysis for PD risk (https://www.nature.com/articles/s41588-023-01584-8) 

<details>
    <summary>MR-MEGA first run logs</summary>
            
```python=
###################
# MR-MEGA v.0.2
###################

Using following command line options:
Input file: MR-MEGA_input.in
Output result file: MAMA_grant.result
Output log file: MAMA_grant.log
Number of PC-s in regression: 5
Binary trait analysis (expecting columns OR, OR_95L, OR_95U)
P-value threshold for showing cohort effect direction: 1
No column filters set
Column names:
	Marker name: MARKERNAME
	Effect allele: EA
	Other allele: NEA
	Effect allele frequency: EAF
	Effect (OR): OR
	Upper CI of effect: OR_95U
	Lower CI of effect: OR_95L
	Strand: STRAND
	Sample size: N
	Chromosome: CHR CHROMOSOME
	Position: POS POSITION
Cohorts list:
saige_phase1.tsv
saige_phase2.tsv
eurGP2_modified.tsv
south_africanGWAS_modified.tsv
africaNo23me_modified.tsv
fooAll_hg38_modified_final.tsv
andrew_modified.tsv
kishore_modified.tsv
Lambda:0.872896
Lambda:0.837092
Lambda:1.3263
Lambda:1.11249
Lambda:1.00554
Lambda:1.09058
Lambda:1.09647
Lambda:1.09521
Principal components:
PCs PC0 PC1 PC2 PC3 PC4
saige_phase1.tsv -0.0411474 0.0337113 -0.0508245 0.0133389 -0.000196743
saige_phase2.tsv -0.0416396 0.032988 -0.0536784 0.0153616 0.000252742
eurGP2_modified.tsv -0.0323588 0.0525199 0.0359604 -0.0476872 -0.00874063
south_africanGWAS_modified.tsv 0.0346749 0.0237841 0.0174668 -0.0220008 0.016018
africaNo23me_modified.tsv 0.155695 -0.00250213 -0.0159387 0.0060501 -0.0056321
fooAll_hg38_modified_final.tsv -0.0306225 -0.10716 -0.0249414 -0.0294943 -0.000376425
adrew_modified.tsv -0.0223473 -0.0171195 0.045521 0.0321814 -0.000507486
kishore_modified.tsv -0.0222544 -0.0162219 0.0464348 0.0322503 -0.000817397

Analysis finished.

```
</details>

### PCA plots 

As it is evidenced in the PCA plots. Together PC1 through PC4 can separate almost all the cohorts, except the South African cohort from the European cohort. These last two are well separated only by PC5. Hence, we must use them all.

Unlike the last multi-ancestry PD meta-analysis where they showed that only 3 PCs were enough to separate all their cohorts, we needed two additional PCs to achieve good separation. However, at that time, neither South Asian or South African cohorts were available. 

<details>
    <summary>PCs</summary>
            
<img width="602" height="490" alt="pc1-2" src="https://github.com/user-attachments/assets/7c380665-e0db-4cf2-9acd-d53e9b017600" />
<img width="602" height="489" alt="pc2-3" src="https://github.com/user-attachments/assets/4cd8f705-0e8c-4670-a898-7b467a033673" />
<img width="602" height="489" alt="pc3-4" src="https://github.com/user-attachments/assets/c1882ef6-8b07-426d-868f-35cc0be5b74e" />
<img width="599" height="489" alt="pc4-5" src="https://github.com/user-attachments/assets/da64147f-c118-4034-8c07-478e024cb7d8" />

</details>


## 3) MR-MEGA second round run

Considering that in the first round we concluded that at least 5 PCs were needed to achieve good separation among the different cohorts, and that these were the defaluts. We did not need a second round run. 

For the record, here is the command line used to run MR-MEGA:
```
MR-MEGA -i MR-MEGA_input.in --pc 5 --o MAMA_grant --name_chr CHR --name_pos POS
```
The MR-MEGA_input.in contains the list with the paths for each of the eight GWASs.

## 4) Format the output summary statistics 

Considering that the raw output from MR-MEGA looks like this:
<details>
    <summary>raw MR-MEGA output header</summary>
            
```python=
MarkerName	Chromosome	Position	EA	NEA	EAF	Nsample	Ncohort	Effects	beta_0	se_0	beta_1	se_1	beta_2	se_2	beta_3	se_3	beta_4	se_4	beta_5	se_5	chisq_association	ndf_association	P-value_association	chisq_ancestry_het	ndf_ancestry_het	P-value_ancestry_het	chisq_residual_het	ndf_residual_het	P-value_residual_het	lnBF	Comments
chr1:730869:C:T	1	730869	T	C	0.0121143	1.84695e+06	5	+-+-???+	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	SmallCohortCount
chr1:758443:G:C	1	758443	C	G	0.117211	1.85059e+06	6	-----??+	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	SmallCohortCount
chr1:762107:A:C	1	762107	C	A	0.0117554	1.83426e+06	3	+++?????

```
</details>

A couple of steps were taken in order to plot the summary statistics, perform the bayesian fine mapping of the genomic risk loci, and the assessment of the discovery of potential novel loci. 

As a first formatting step, we use some help from the python toolkit called gwaslab (https://cloufield.github.io/gwaslab/). 

We use their function called basic_check() to get an initial clean summary statistics from which we could further elaborate on.

<details>
    <summary>first_pass_parse_gwaslab.py</summary>
            
```python=
import gwaslab as gl

input_file="MAMA_grant.result"
output_file="MR-MEGA"
cores=12

ss = gl.Sumstats(input_file, build="38",
    snpid="MarkerName",     
    chrom="Chromosome",
    pos="Position",
    ea="EA",
    eaf="EAF",
    nea="NEA",
    n="Nsample",
    p="P-value_association",
    direction="Effects",
    phet="P-value_ancestry_het",
    other=["Ncohort","P-value_residual_het","lnBF"],
    sep="\t",
    na_values=["NA","."],
    verbose=True
)

ss.basic_check(
    remove=True,
    fixid_args={"fixchrpos":False,"fixid":False,"fixsep":False,"forcefixid":False,"overwrite":False},
    n_cores=cores,
    removedup_args={"mode":"c","keep":"first","keep_col":"P","remove":True}, 
    normalizeallele_args={},
    verbose=True
)

ss.to_format(path=output_file, fmt="gwaslab",cols=["Ncohort","P-value_residual_het","lnBF"])

```
</details>

Once this file is generated, we call gwaslab again to assign each entry an rsID based on the most recent dbSNP reference files for build hg38. 
This functionality could have been included in the previous step, however it is nice to have a clean copy of your data. Gwaslab offers as well some functions to save intermediate summary statistics files. This could be another alternative. 
 
<details>
    <summary>assign_rsID.py</summary>
            
```python=
import gwaslab as gl

input_file="MR-MEGA.gwaslab.tsv" #the file that the first parse outputs 
output_file="MR-MEGA_with_rsID_fullCols_hg38"
logFile="getrsID_log"
cores=12

#Load the data 
ss = gl.Sumstats(input_file, build="38", #you could potentially input the data as gwaslab format, however, if you do not specify which additional columns to be included here and at the output, they are going to be ignored
    snpid="SNPID",     
    chrom="CHR",
    pos="POS",
    ea="EA",
    eaf="EAF",
    nea="NEA",
    n="N",
    p="P",
    direction="DIRECTION",
    phet="P_HET",
    other=["Ncohort","P-value_residual_het","lnBF"],
    sep="\t",
    na_values=["NA","."],
    verbose=True
)

ss.basic_check(
    remove=False,
    n_cores=cores,
    normalizeallele_args={},
    verbose=True
)

ss.assign_rsid(n_cores=cores,
                       ref_rsid_vcf="/home/duartej3/beegfs/JF/programs/gwaslab/referenceFiles/GCF_000001405.40.gz",
                       chr_dict = gl.get_number_to_NC(build="38") 
)


ss.sort_column()
ss.to_format(path=output_file, fmt="gwaslab",cols=["Ncohort","P-value_residual_het","lnBF"])
ss.log.save(path=logFile)

```
</details>

The output summary statistics with the rsIDs look like this (P_HET is the column that contains the p value for the ancestry heterogeneity test):

```
SNPID	rsID	CHR	POS	EA	NEA	EAF	P	N	DIRECTION	P_HET	STATUS	Ncohort	P-value_residual_het	lnBF
chr1:858952:G:A	rs12127425	1	858952	A	G	0.0771	9.1735e-01	1884050	-----+++	0.944081	3860099	8	0.410831	-5.22565
chr1:910255:C:T	rs117086422	1	910255	T	C	0.197	2.2446e-01	1884050	++--++--	0.151687	3860099	8	0.803276	-2.14298
```

## 5) Fine mapping and novel genomic risk loci 

To fine map the genome-wide significant genomic risk loci we employed the same framework used in the most recent multi-ancestry meta-analysis for PD risk (https://www.nature.com/articles/s41588-023-01584-8) and our LARGE-PD gwas (https://www.medrxiv.org/content/10.1101/2025.07.18.25331793v1.full-text) That leverages the Bayesian factor outputted by MR-MEGA to construct 95 and 99% credible sets for each genomic region discovered. 

This requires uploading the summary statistics to the FUMA platform (https://fuma.ctglab.nl/). For the upload settings, we choose a p-value threshold of 5e-9 to define the significant genomic risk loci. Since the increased number of ancestries in the dataset potentially increases the number of haplotypes/number of independent tests. For the LD clumping to define the lead and independent variants, we choose the 1KGP “all” reference panel. 

Since FUMA requires either the rsIDs or the genomic positions to be in hg19, we performed a liftover of the sumstats to hg19. Specially since the fine mapping pipeline relies on the genomic risk loci defined by FUMA in order to perform bayesian fine mapping. 

Briefly, the script loads the FUMA GenomicRiskLoci.txt file, which defines each genomic risk locus by chromosome and hg19 start/end positions. For each risk locus, we subset the MR-MEGA summary statistics to variants within that interval and convert MR-MEGA’s log Bayes factors (lnBF) to Bayes factors (BF = exp(lnBF)) for association versus the null model. Assuming a single causal variant per locus and equal prior probability across SNPs, we treat the Bayes factors as proportional to the posterior probability of causality, normalize them to obtain posterior probabilities, and rank SNPs accordingly. We then build 95% and 99% credible sets by adding SNPs in order of decreasing posterior probability until the cumulative posterior mass within the locus reaches 95% or 99%, respectively. These credible sets and their posterior probabilities are recorded in the fine-mapping report.

<details>
    <summary>bayesian_fine_mapping.py</summary>
            
```python=
import os
import argparse
import pandas as pd
import numpy as np
from scipy.stats import chi2


def load_fuma_and_summary(fuma_dir: str, mrmega_file: str):
    #Load FUMA risk-locus definitions and MR-MEGA summary statistics.
    loci_path = os.path.join(fuma_dir, "GenomicRiskLoci.txt")
    if not os.path.exists(loci_path):
        raise FileNotFoundError(f"Cannot find FUMA loci file at {loci_path}")
    loci = pd.read_csv(loci_path, delim_whitespace=True)

    if not os.path.exists(mrmega_file):
        raise FileNotFoundError(f"Cannot find MR-MEGA file at {mrmega_file}")
    sumstat = pd.read_csv(
        mrmega_file,
        delim_whitespace=True,
        compression='infer',
        engine='c'
    )
    return loci, sumstat


def finemap_locus(loci_df, sumstat_df):
    """
    Compute 95% and 99% credible sets for each locus based on Bayes Factors.
    Returns report, cs95_df, cs99_df.
    """

        # Work on copies so we don't mutate original dataframes outside
    loci_df = loci_df.copy()
    sumstat_df = sumstat_df.copy()

    # Harmonize chromosome types 
    if 'CHR' in sumstat_df.columns:
        sumstat_df['CHR'] = pd.to_numeric(sumstat_df['CHR'], errors='coerce')
    if 'chr' in loci_df.columns:
        loci_df['chr'] = pd.to_numeric(loci_df['chr'], errors='coerce')

    # Harmonize positions 
    for col in ['start', 'end']:
        if col in loci_df.columns:
            loci_df[col] = pd.to_numeric(loci_df[col], errors='coerce')

    if 'POS_HG37' not in sumstat_df.columns:
        raise KeyError("Expected column 'POS_HG37' in MR-MEGA sumstats.")

    sumstat_df['POS_HG37'] = pd.to_numeric(sumstat_df['POS_HG37'], errors='coerce')
    
    # Convert lnBF to BF
    sumstat_df['BF'] = np.exp(sumstat_df['lnBF'])
    report_rows = []
    cs95_list = []
    cs99_list = []
    grouped = {chrom: df for chrom, df in sumstat_df.groupby('CHR')}

    for idx, row in loci_df.iterrows():
        chrom = row['chr']
        start = row['start']
        end = row['end']
        df_chr = grouped.get(chrom)
        if df_chr is None:
            continue
        window = df_chr[(df_chr['POS_HG37'] >= start) & (df_chr['POS_HG37'] <= end)].copy()
        window = window.sort_values('BF', ascending=False).reset_index(drop=True)
        total_bf = window['BF'].sum()

        # 95% credible set
        cum = 0.0
        sel95 = []
        for i, bf in enumerate(window['BF']):
            cum += bf
            sel95.append(i)
            if cum / total_bf >= 0.95:
                break
        cs95 = window.loc[sel95].copy()
        cs95['PP'] = cs95['BF'] / total_bf
        cs95['Locus'] = idx + 1
        cs95_list.append(cs95)

        # 99% credible set
        cum = 0.0
        sel99 = []
        for i, bf in enumerate(window['BF']):
            cum += bf
            sel99.append(i)
            if cum / total_bf >= 0.99:
                break
        cs99 = window.loc[sel99].copy()
        cs99['PP'] = cs99['BF'] / total_bf
        cs99['Locus'] = idx + 1
        cs99_list.append(cs99)

        # Collect SNP lists per locus
        snps95 = window.loc[sel95, 'MarkerName'] if 'MarkerName' in window.columns else window.loc[sel95, 'rsID']
        snps99 = window.loc[sel99, 'MarkerName'] if 'MarkerName' in window.columns else window.loc[sel99, 'rsID']

        report_rows.append({
            'Locus': idx + 1,
            'NumSNPs': window.shape[0],
            'NumSNPs_in_95%': len(sel95),
            'NumSNPs_in_99%': len(sel99),
            'SNPs_in_95%': ';'.join(snps95.astype(str)),
            'SNPs_in_99%': ';'.join(snps99.astype(str))
        })

    report = pd.DataFrame(report_rows)
    cs95_df = pd.concat(cs95_list, ignore_index=True) if cs95_list else pd.DataFrame()
    cs99_df = pd.concat(cs99_list, ignore_index=True) if cs99_list else pd.DataFrame()
    return report, cs95_df, cs99_df


def annotate_with_fuma(cs_df, fuma_dir: str):
    """
    Annotate a credible-set DataFrame with FUMA snp annotations.
    """
    snp_file = os.path.join(fuma_dir, 'snps.txt')
    if not os.path.exists(snp_file):
        raise FileNotFoundError(f"Cannot find FUMA SNP file at {snp_file}")
    fuma = pd.read_csv(snp_file, delim_whitespace=True)
    fuma = fuma.rename(columns={'chr': 'CHR', 'pos': 'BP'})
    annotated = cs_df.merge(
        fuma,
        left_on=['CHR', 'POS_HG37'],
        right_on=['CHR', 'BP'],
        how='left'
    )
    return annotated


def main():
    parser = argparse.ArgumentParser(description="Fine-map MR-MEGA results using FUMA loci.")
    parser.add_argument('--fuma-dir', required=True, help='Path to FUMA output directory')
    parser.add_argument('--mrmega-file', required=True, help='Path to MR-MEGA summary stats (.tsv.gz)')
    parser.add_argument('--out-dir', default='output', help='Directory for fine-mapping outputs')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    loci, sumstat = load_fuma_and_summary(args.fuma_dir, args.mrmega_file)
    report, cs95, cs99 = finemap_locus(loci, sumstat)

    # Save reports
    report_file = os.path.join(args.out_dir, 'finemap_report.tsv')
    cs95_file   = os.path.join(args.out_dir, 'finemap_cs_95.tsv')
    cs99_file   = os.path.join(args.out_dir, 'finemap_cs_99.tsv')
    report.to_csv(report_file, sep='\t', index=False)
    cs95.to_csv(cs95_file, sep='\t', index=False)
    cs99.to_csv(cs99_file, sep='\t', index=False)

    # Annotate credible sets
    annotated95 = annotate_with_fuma(cs95, args.fuma_dir)
    annotated99 = annotate_with_fuma(cs99, args.fuma_dir)
    annotated95_file = os.path.join(args.out_dir, 'finemap_cs_95_annot.tsv')
    annotated99_file = os.path.join(args.out_dir, 'finemap_cs_99_annot.tsv')
    annotated95.to_csv(annotated95_file, sep='\t', index=False)
    annotated99.to_csv(annotated99_file, sep='\t', index=False)

    # Identify singleton loci and their rsIDs
    single95_loci = report[report['NumSNPs_in_95%'] == 1]['Locus']
    single95_snps = cs95[cs95['Locus'].isin(single95_loci)][['Locus', 'rsID']]
    single95_snps = single95_snps.rename(columns={'rsID': 'rsID'})
    single95_file = os.path.join(args.out_dir, 'singleton_loci_95.tsv')
    single95_snps.to_csv(single95_file, sep='\t', index=False)

    single99_loci = report[report['NumSNPs_in_99%'] == 1]['Locus']
    single99_snps = cs99[cs99['Locus'].isin(single99_loci)][['Locus', 'rsID']]
    single99_snps = single99_snps.rename(columns={'rsID': 'rsID'})
    single99_file = os.path.join(args.out_dir, 'singleton_loci_99.tsv')
    single99_snps.to_csv(single99_file, sep='\t', index=False)

if __name__ == "__main__":
    main()

```
</details>

An example on how to run the script is as follow:

```
python fine_map.py --fuma-dir ./ --mrmega-file MR-MEGA_with_rsID_fullCols_hg37.tsv
```

The fine mapping results could be found in the data repository accompaying this github. 
We report 11 genomic risk regions that have a single variant in the 99% credible set and 13 genomic risk regions that have a single variant in the 95% credible set. 


## 6) Interrogation of novel risk loci 

In order to interrogate if we have genome-wide significant novel loci, we have to gather all the previously known loci associated with Parkinson’s Disease. 

For this purpose, we downloaded all the associations reported by GWAS catalog (https://www.ebi.ac.uk/gwas/efotraits/MONDO_0005180) and the NDKP database (https://ndkp.hugeamp.org/phenotype.html
phenotype=Parkinsons). 

GWAS catalog stores all the variants associated with PD with a p value of at least 1e-5 within their datasets. While NDKP reports all the associated variants with a p value of at least 5e-8 within their
datasets. Hence, in order to be as comprehensive as possible, we gathered all the rsID available in both databases (over 800 in GWAS catalog and over 200 in NDKP) and made a single text list (attached with the
data repository accompaying the github). 

With this rsID table, we leveraged a gwaslab lab function to assign chromosome and position in build hg38 to each of the SNPs. We did it this way in order to avoid potential discordances among the genomic
coordinates that GWAS catalog and NDKP reported, since there was ambiguity on the build of the positions. Hence, using their own rsIDs to map the adequate genomic coordinates is the most reliable procedure. 

The chromosome and position of each rsID will be used by gwaslab again to compare against our MR-MEGA summary statistics and define novel risk loci. 

With the following code we obtained the CHR and POS for each rsID using a modified version of the dbSNP reference file (build hg38)

<details>
    <summary>get_rsIDs_chr_pos.py</summary>
            
```python=
import gwaslab as gl

input_file="gwasCatalog_ndkp_known_rsID_PD.txt"
output_file="gwasCatalog_ndkp_known_rsID_PD_chr_pos"
logFile="gwasCatalog_ndkp_known_rsID_PD_chr_pos"
cores=12

ss = gl.Sumstats(input_file, build="38",
    rsid="rsID", 
    verbose=True
)

ss.rsid_to_chrpos2(path="/home/duartej3/beegfs/JF/programs/gwaslab/referenceFiles/GCF_000001405.40.gz.rsID_CHR_POS_groups_20000000.h5",
			#n_cores="12",
			build="38"
)


ss.fix_chr()
ss.fix_pos()
ss.sort_coordinate()
ss.sort_column()
ss.to_format(path=output_file, fmt="gwaslab")
ss.log.save(path=logFile)
```
</details>

The file obtained looks like this:

```
rsID	CHR	POS	STATUS
rs7532024	1	7119419	3890999
rs302714	1	8426071	3890999
```

The official table that gwaslab uses to compare against in order to find novel loci must have only chromosome and position and must be free of NA values, like this:
```
CHR	POS
1	7119419
1	8426071
```

Using gwaslab, we leveraged the known loci table that we just built and interrogated the presence of novel risk loci. Based on the criteria outlined in the latest multi-ancestry meta-analysis for PD risk (https://www.nature.com/articles/s41588-023-01584-8) 

<details>
    <summary>get_novelLoci.py</summary>
            
```python=
import gwaslab as gl

input_file="MR-MEGA_with_rsID_fullCols_hg38"
logFile="novel_log"
cores=12

ss = gl.Sumstats(input_file, build="38", fmt="gwaslab",

)

#new
ss.fix_chr(remove=True)
ss.fix_pos(remove=True)

res = ss.get_novel(
    known="gwasCatalog_ndkp_known_rsID_PD_just_chr_pos.tsv",
    only_novel=False,            # return both
    sig_level=5e-9,
    windowsizekb=250,
    windowsizekb_for_novel=250
)
res.query("NOVEL == True").to_csv("novel_hits.tsv", sep="\t", index=False)
res.query("NOVEL == False").to_csv("known_hits.tsv", sep="\t", index=False)

ss.log.save(path=logFile)
```
</details>

Here is the log file, in which is registered the discovery of 8 potential novel loci:

<details>
    <summary>log file novel loci</summary>
            
```python=
2025/11/20 21:16:17 Sumstats Object created.
2025/11/20 21:16:17 GWASLab v3.6.10 https://cloufield.github.io/gwaslab/
2025/11/20 21:16:17 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/11/20 21:16:17 Python version: 3.12.2 | packaged by conda-forge | (main, Feb 16 2024, 20:50:58) [GCC 12.3.0]
2025/11/20 21:16:17 Start to load format from formatbook....
2025/11/20 21:16:17  -gwaslab format meta info:
2025/11/20 21:16:17   - format_name  : gwaslab
2025/11/20 21:16:17   - format_source  : https://cloufield.github.io/gwaslab/
2025/11/20 21:16:17   - format_version  : 20231220_v4
2025/11/20 21:16:17 Start to initialize gl.Sumstats from file :../liftoverMR-MEGA_rsID_noSNPID_sumstats/MR-MEGA_with_rsID_fullCols_noSNPID.gwaslab.tsv
2025/11/20 21:16:20  -Reading columns          : EA,P,N,EAF,P_HET,rsID,CHR,STATUS,POS,NEA,DIRECTION
2025/11/20 21:16:20  -Renaming columns to      : EA,P,N,EAF,P_HET,rsID,CHR,STATUS,POS,NEA,DIRECTION
2025/11/20 21:16:20  -Current Dataframe shape : 2278131  x  11
2025/11/20 21:16:20  -Initiating a status column: STATUS ...
2025/11/20 21:16:20  -Genomic coordinates are based on GRCh38/hg38...
2025/11/20 21:16:20 Start to reorder the columns...v3.6.10
2025/11/20 21:16:20  -Current Dataframe shape : 2278131 x 11 ; Memory usage: 174.31 MB
2025/11/20 21:16:20  -Reordering columns to    : rsID,CHR,POS,EA,NEA,EAF,P,N,DIRECTION,P_HET,STATUS
2025/11/20 21:16:20 Finished reordering the columns.
2025/11/20 21:16:20  -Trying to convert datatype for CHR: string -> Int64...Int64
2025/11/20 21:16:21  -Column  : rsID   CHR   POS   EA       NEA      EAF     P       N     DIRECTION P_HET   STATUS  
2025/11/20 21:16:21  -DType   : object Int64 int64 category category float64 float64 int64 object    float64 category
2025/11/20 21:16:21  -Verified: T      T     T     T        T        T       T       T     T         T       T       
2025/11/20 21:16:21  -Current Dataframe memory usage: 176.48 MB
2025/11/20 21:16:21 Finished loading data successfully!
2025/11/20 21:16:21  -Genomic coordinates are based on GRCh38/hg38...
2025/11/20 21:16:21 Path component detected: ['23453636628688']
2025/11/20 21:16:21 Creating path: ./23453636628688
2025/11/20 21:16:21 Start to fix chromosome notation (CHR)...v3.6.10
2025/11/20 21:16:21  -Current Dataframe shape : 2278131 x 11 ; Memory usage: 176.48 MB
2025/11/20 21:16:21  -Checking CHR data type...
2025/11/20 21:16:21  -Variants with standardized chromosome notation: 2278131
2025/11/20 21:16:21  -All CHR are already fixed...
2025/11/20 21:16:23 Finished fixing chromosome notation (CHR).
2025/11/20 21:16:23 Start to fix basepair positions (POS)...v3.6.10
2025/11/20 21:16:23  -Current Dataframe shape : 2278131 x 11 ; Memory usage: 176.48 MB
2025/11/20 21:16:23  -Converting to Int64 data type ...
2025/11/20 21:16:25  -Position bound:(0 , 250,000,000)
2025/11/20 21:16:25  -Removed outliers: 0
2025/11/20 21:16:25  -Removed 0 variants with bad positions.
2025/11/20 21:16:25 Finished fixing basepair positions (POS).
2025/11/20 21:16:25 Start to check if lead variants are known...v3.6.10
2025/11/20 21:16:25  -Current Dataframe shape : 2278131 x 11 ; Memory usage: 196.04 MB
2025/11/20 21:16:25 Start to extract lead variants...v3.6.10
2025/11/20 21:16:25  -Current Dataframe shape : 2278131 x 11 ; Memory usage: 196.04 MB
2025/11/20 21:16:25  -Processing 2278131 variants...
2025/11/20 21:16:25  -Significance threshold : 5e-09
2025/11/20 21:16:25  -Sliding window size: 250  kb
2025/11/20 21:16:26  -Using P for extracting lead variants...
2025/11/20 21:16:26  -Found 4123 significant variants in total...
2025/11/20 21:16:26  -Identified 70 lead variants!
2025/11/20 21:16:26 Finished extracting lead variants.
2025/11/20 21:16:26  -Lead variants in known loci: 1033
2025/11/20 21:16:26  -Checking the minimum distance between identified lead variants and provided known variants...
2025/11/20 21:16:26  -Identified  62  known vairants in current sumstats...
2025/11/20 21:16:26  -Identified  8  novel vairants in current sumstats...
2025/11/20 21:16:26 Finished checking if lead variants are known.
```
</details>

The summary statistics for the known and novel variants identified are in the data repository accompanying the github.

## 7) Plotting 

### Novel loci 
<img width="7466" height="2766" alt="plot_novelLoci" src="https://github.com/user-attachments/assets/15c180d2-6b61-4aa2-92cc-37f7ebf7d7f5" />

We used this script to highlight the 8 novel loci (in red and with their respective nearest protein coding gene):
<details>
    <summary>plot_novelLoci.py</summary>
            
```python=
import gwaslab as gl

input_file="MR-MEGA_with_rsID_fullCols_hg38"
outPlot="plot_novelLoci.png"
logFile="plot_novelLoci.txt"

ss = gl.Sumstats(input_file, build="38", fmt="gwaslab"
)


ss.plot_mqq(
                  mode="m",
                  sig_level=5e-9,
                  build="38", 
                  anno="GENENAME",
                  anno_set=["rs269291","rs4954490","rs2445964","rs6578471","rs933738","rs75544157","rs2837732","rs74793276"] ,
                  pinpoint=["rs269291","rs4954490","rs2445964","rs6578471","rs933738","rs75544157","rs2837732","rs74793276"],
                  stratified=True,
                  marker_size=(5,10),
                  save=outPlot, save_args={"dpi":600,"facecolor":"white"}
                  )

ss.log.save(path=logFile)
```
</details>

### Overall plot 

<img width="7430" height="2716" alt="plot" src="https://github.com/user-attachments/assets/b902b62d-cf17-4c97-85f0-c30af1bf8171" />

We used this script to plot the 70 loci that were genome-wide significant at the multi-ancestry level (5-e9).

<details>
    <summary>plot_overallLoci.py</summary>
            
```python=
import gwaslab as gl

input_file="MR-MEGA_with_rsID_fullCols_hg38"
outPlot="plot_novelLoci.png"
logFile="plot_novelLoci.txt"

ss = gl.Sumstats(input_file, build="38", fmt="gwaslab"
)


ss.plot_mqq(
                  mode="m",
                  sig_level=5e-9,
                  repel_force=0.1,
                  build="38", 
                  anno="GENENAME",
                  anno_style="tight",
                  anno_fontsize=5,
                  font_family="DejaVu Sans",
                  #anno_set=["rs269291","rs4954490","rs2445964","rs6578471","rs933738","rs75544157","rs2837732","rs74793276"] ,
                  #pinpoint=["rs269291","rs4954490","rs2445964","rs6578471","rs933738","rs75544157","rs2837732","rs74793276"],
                  stratified=True,
                  marker_size=(5,10),
                  save=outPlot, save_args={"dpi":600,"facecolor":"white"}
                  )

ss.log.save(path=logFile)
```
</details>




