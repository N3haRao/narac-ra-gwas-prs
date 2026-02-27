# Load required libraries
library(GMMAT)

# Load data
pheno <- read.table("/projectnb/bs859/students/neharao/final-project/ra_cleaned.fam", header = FALSE)
colnames(pheno) <- c("FID", "IID", "fa", "mo", "sex", "case")

pcs <- read.table("/projectnb/bs859/students/neharao/final-project/RA_pcs.txt", header = TRUE)

# Filter females only
pheno_female <- pheno[pheno$sex == 2, ]  # assuming 2 represents females

# Merge PC data with fam file data
pheno_female_pcs <- merge(pheno_female, pcs, by.x = c("FID", "IID"), by.y = c("FID", "IID"), all.x = TRUE)

# Load GRM (genetic relationship matrix)
grm <- as.matrix(read.table("q2a-grm-female.rel", header = FALSE))
grm_ids <- read.table("q2a-grm-female.rel.id", header = FALSE)

# Apply IDs to GRM matrix
dimnames(grm)[[1]] <- dimnames(grm)[[2]] <- grm_ids[, 2]

# Define models
female_model_null <- glmmkin(case - 1 ~ 1, data = pheno_female_pcs, id = "IID", kins = grm, family = binomial("logit"))
female_model_covariates <- glmmkin(case - 1 ~ PC2 + PC3, data = pheno_female_pcs, id = "IID", kins = grm, family = binomial("logit"))

# Perform genome-wide association analysis
geno.file <- "/projectnb/bs859/students/neharao/final-project/ra_cleaned"
glmm.score(female_model_null, infile = geno.file, outfile = "test.glmm.score.female.null")
glmm.score(female_model_covariates, infile = geno.file, outfile = "test.glmm.score.female.covariates")
