# 🧬 narac-ra-gwas-prs

### Sex-Stratified GWAS, Population Stratification Control, and Polygenic Risk Scoring in Rheumatoid Arthritis (NARAC)

**Author:** Neha Rao

**Tools:** PLINK · EIGENSOFT · GMMAT · PRSice2 · LDSC · R · Bash · Python

**Environment:** HPC (Linux cluster)

**Study Type:** Case-control GWAS with mixed model association testing

---

## 🔎 Project Overview

This repository contains a fully reproducible genome-wide association analysis (GWAS) pipeline for **Rheumatoid Arthritis (RA)** using data from the **North American Rheumatoid Arthritis Consortium (NARAC)**.

The analysis includes:

* Rigorous genotype quality control
* Population structure correction via PCA
* Sex-stratified mixed-model GWAS
* Polygenic Risk Score (PRS) analysis using external summary statistics
* SNP heritability estimation using LD Score Regression

This project demonstrates end-to-end genetic epidemiology workflow design, statistical modeling, and HPC-based genomic data processing.

---

## 🧬 Dataset

**North American Rheumatoid Arthritis Consortium (NARAC)**

* 868 RA cases
* 1,194 controls
* Predominantly Northern European ancestry
* Case-control design
* Unrelated individuals

Cases met American College of Rheumatology criteria. Controls were sourced from the New York Cancer Project.

**Note:** NARAC genotype data are controlled-access and not included in this repository. This repo provides code and workflow only.

Acknowledgment: Genetic Analysis Workshop 16 (GAW16)
[https://doi.org/10.1002/gepi.20464](https://doi.org/10.1002/gepi.20464)

---

# 📊 Analysis Workflow

Below is the full analytical pipeline implemented in this project.

---

## 1️⃣ Genotype Quality Control (PLINK)

**SNP-Level Filters**

* Minor Allele Frequency (MAF) > 0.01
* Missingness (geno) < 5%
* Hardy–Weinberg Equilibrium p > 1e-6

**Sample-Level Filters**

* Individual missingness (mind) < 5%

Output:

```
cleaned/ra_cleaned.*
```

---

## 2️⃣ LD Pruning + PCA (EIGENSOFT smartpca)

To control for population stratification:

* LD pruning (`--indep-pairwise`)
* PCA on pruned SNP set
* Extraction of top 10 PCs
* Logistic regression of PCs vs case status to identify significant covariates

Output:

```
pca/ra_pruned.evec
pca/RA_pcs.txt
```

---

## 3️⃣ Sex-Stratified GWAS (Mixed Model — GMMAT)

To account for relatedness and residual structure:

* Built Genetic Relationship Matrix (GRM) per sex
* Fitted logistic mixed models (`glmmkin`)
* Genome-wide score test (`glmm.score`)

### Models:

**Null model**

```
case ~ 1
```

**Covariate model**

```
case ~ PC2 + PC3
```

(PCs selected based on association testing)

Output:

```
gwas/test.glmm.score.female.*
gwas/test.glmm.score.male.*
```

---

## 4️⃣ Polygenic Risk Scores (PRSice2)

PRS computed using external GWAS summary statistics:

* European meta-analysis
* Asian meta-analysis

Covariates included selected PCs.

Output:

```
prs/PRS_European.*
prs/PRS_Asian.*
```

---

## 5️⃣ SNP Heritability (LD Score Regression)

* Summary statistics munged using `munge_sumstats.py`
* Heritability estimated using ancestry-matched LD scores
* Liability-scale correction applied

Output:

```
RA_h2_European.log
RA_h2_Asian.log
```

---

# 📈 Key Methodological Strengths

✔ Mixed-model GWAS accounting for relatedness
✔ Explicit correction for population stratification
✔ Sex-stratified analysis
✔ Cross-ancestry PRS comparison
✔ HPC reproducible pipeline design
✔ Integration of Bash, R, Python, and external genomic tools

---

# 🧠 Technical Skills Demonstrated

* Genome-wide association analysis
* Logistic mixed modeling
* Genetic relationship matrices
* Polygenic risk scoring
* LD score regression
* Linux/HPC scripting
* End-to-end genomic pipeline design
* R statistical modeling
* Python file parsing
* Reproducible research practices

---

# 📚 References

- Purcell S, Neale B, Todd-Brown K, et al. PLINK: A tool set for whole-genome association and population-based linkage analyses. *Am J Hum Genet*. 2007;81(3):559-575. [doi:10.1086/519795](https://doi.org/10.1086/519795)
- Price AL, Patterson NJ, Plenge RM, et al. Principal components analysis corrects for stratification in genome-wide association studies. *Nat Genet*. 2006;38(8):904-909. [doi:10.1038/ng1847](https://doi.org/10.1038/ng1847)
- Anderson CA, Pettersson FH, Clarke GM, et al. Data quality control in genetic case-control association studies. *Nat Protoc*. 2010;5(9):1564-1573. [doi:10.1038/nprot.2010.116](https://doi.org/10.1038/nprot.2010.116)
- Wray NR, Goddard ME, Visscher PM. Prediction of individual genetic risk to disease from genome-wide association studies. *Genome Res*. 2007;17(10):1520-1528. [doi:10.1101/gr.6665407](https://doi.org/10.1101/gr.6665407)

