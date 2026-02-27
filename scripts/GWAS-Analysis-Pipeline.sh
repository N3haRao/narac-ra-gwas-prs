#!/bin/bash

# Final Project Code: NARAC GWAS Analysis
# Author: Neha Rao
# Purpose: Conduct GWAS analysis for the NARAC dataset, including data cleaning, PCA, GWAS, and PRS.

# Define base paths
BASE_DIR="/projectnb/bs859/students/neharao/final-project"
CLEANED_DIR="${BASE_DIR}/cleaned"
PRUNED_DIR="${BASE_DIR}/pruned"
PCA_DIR="${BASE_DIR}/pca"
GWAS_DIR="${BASE_DIR}/gwas"
PRS_DIR="${BASE_DIR}/prs"

mkdir -p $CLEANED_DIR $PRUNED_DIR $PCA_DIR $GWAS_DIR $PRS_DIR

# Section 1: Genetic Data Cleaning
# Step 1: Filter SNPs based on minor allele frequency (MAF), genotype call rate, and Hardy-Weinberg equilibrium (HWE)
echo "Step 1: SNP Quality Control"
plink --bfile ${BASE_DIR}/RA-fp/narac_hg19 \
  --maf 0.01 --geno 0.05 --hwe 1e-6 --freq --make-bed --out ${CLEANED_DIR}/ra_cleaned_1

# Step 2: Remove individuals with high missing genotype rates
echo "Step 2: Individual Quality Control"
plink --bfile ${CLEANED_DIR}/ra_cleaned_1 \
  --mind 0.05 --make-bed --out ${CLEANED_DIR}/ra_cleaned

# Section 2: Pruning and PCA
# Step 1: Prune markers for PCA (removes LD-correlated SNPs)
echo "Step 3: Pruning SNPs"
plink --bfile ${CLEANED_DIR}/ra_cleaned \
  --geno 0.01 --maf 0.02 --indep-pairwise 10000kb 1 0.15 --out ${PRUNED_DIR}/ra

# Step 2: Apply pruning results
echo "Step 4: Generating Pruned Dataset"
plink --bfile ${CLEANED_DIR}/ra_cleaned \
  --extract ${PRUNED_DIR}/ra.prune.in --make-bed --out ${PRUNED_DIR}/ra_cleaned_pruned

# Step 3: Perform PCA
echo "Step 5: Running PCA"
cat > ${PCA_DIR}/q1.par <<EOF
genotypename: ${PRUNED_DIR}/ra_cleaned_pruned.bed
snpname: ${PRUNED_DIR}/ra_cleaned_pruned.bim
indivname: ${PRUNED_DIR}/ra_cleaned_pruned.fam
evecoutname: ${PCA_DIR}/ra_pruned.evec
evaloutname: ${PCA_DIR}/ra_pruned.eval
altnormstyle: NO
numoutevec: 10
numoutlieriter: 0
outliersigmathresh: 4
outlieroutname: ${PCA_DIR}/outliers.removed
EOF
smartpca -p ${PCA_DIR}/q1.par > ${PCA_DIR}/q1.out

#Step 4: Usage of the file RA_pcs.txt (RA_pcs.txt is the output from smartpca with a column header added) to determine which PCs are associated with case status.
#“--logistic no-snp” produces results for the specified covariates, with no SNPs in the model.
echo "Step 6: Determining PCs associated with case status"
plink --bfile ${CLEANED_DIR}/ra_cleaned --covar ${PCA_DIR}/RA_pcs.txt \
  --covar-name PC1-PC10 --out checkPCs --logistic no-snp beta

# Section 3: GWAS Analysis
# Female GWAS Analysis
echo "Step 7: Running Female GWAS"
plink --bfile ${CLEANED_DIR}/ra_cleaned \
  --filter-females --make-rel square --out ${GWAS_DIR}/q2a-grm-female
#awk commands to replace the "PVAL" column header with "P", "SCORE" with "BETA" and the "POS" with "BP" to match PLINK output
awk 'NR==1{$4="BP";$9="BETA";$11="P"};{print $0}' test.glmm.score.female.nocov>test.glmm.score.female.nocov.txt
awk 'NR==1{$4="BP";$9="BETA";$11="P"};{print $0}' test.glmm.score.female.covariates>test.glmm.score.female.covariates.txt

# Male GWAS Analysis
echo "Step 8: Running Male GWAS"
plink --bfile ${CLEANED_DIR}/ra_cleaned \
  --filter-males --make-rel square --out ${GWAS_DIR}/q2a-grm-male
#awk commands to replace the "PVAL" column header with "P", "SCORE" with "BETA" and the "POS" with "BP" to match PLINK output
awk 'NR==1{$4="BP";$9="BETA";$11="P"};{print $0}' test.glmm.score.male.nocov>test.glmm.score.male.nocov.txt
awk 'NR==1{$4="BP";$9="BETA";$11="P"};{print $0}' test.glmm.score.male.covariates>test.glmm.score.male.covariates.txt

# Section 4: Polygenic Risk Scores (PRS)
# PRS for European Population
echo "Step 9: Generating PRS for European Population"
Rscript $SCC_PRSICE_BIN/PRSice.R --prsice $SCC_PRSICE_BIN/PRSice \
  --base ${BASE_DIR}/RA_GWASmeta_European_v2.txt \
  --target ${CLEANED_DIR}/ra_cleaned \
  --snp SNPID --A1 A1 --A2 A2 --bp BP --stat OR --pvalue P-val \
  --cov-file ${PCA_DIR}/RA_pcs.txt --cov-col PC1,PC4 \
  --binary-target T --out ${PRS_DIR}/PRS_European

# PRS for Asian Population
echo "Step 9: Generating PRS for Asian Population"
Rscript $SCC_PRSICE_BIN/PRSice.R --prsice $SCC_PRSICE_BIN/PRSice \
  --base ${BASE_DIR}/RA_GWASmeta_Asian_v2.txt \
  --target ${CLEANED_DIR}/ra_cleaned \
  --snp SNPID --A1 A1 --A2 A2 --bp BP --stat OR --pvalue P-val \
  --cov-file ${PCA_DIR}/RA_pcs.txt --cov-col PC1,PC4 \
  --binary-target T --out ${PRS_DIR}/PRS_Asian

echo "Analysis Complete. Check output files in the respective directories."
