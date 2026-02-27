#!/bin/bash
module load ldsc
# Set paths to data and LD score files
RA_DATA_DIR="/projectnb/bs859/data/RheumatoidArthritis/final_project"
LDSCORES_DIR="/projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore"

# Prepare summary statistics for European population
munge_sumstats.py --sumstats $RA_DATA_DIR/RA_GWASmeta_European_v2.txt --snp snpid --a1 effect_allele --a2 other_allele --signed-sumstats beta,0 --merge-alleles /projectnb/bs859/data/ldscore_files/w_hm3.snplist --out RA_European

# Prepare summary statistics for Asian population
munge_sumstats.py --sumstats $RA_DATA_DIR/RA_GWASmeta_Asian_v2.txt.gz --snp snpid --a1 effect_allele --a2 other_allele --signed-sumstats beta,0 --merge-alleles /projectnb/bs859/data/ldscore_files/w_hm3.snplist --out RA_Asian

# Prepare summary statistics for transethnic meta-analysis
#munge_sumstats.py --sumstats $RA_DATA_DIR/RA_GWASmeta_TransEthnic_v2.txt.gz --snp snpid --a1 effect_allele --a2 other_allele --signed-sumstats beta,0 --merge-alleles /projectnb/bs859/data/ldscore_files/w_hm3.snplist --out RA_TransEthnic

# Perform LD score regression for European population
ldsc.py --h2 RA_European.sumstats.gz --ref-ld-chr $LDSCORES_DIR/UKBB.EUR.rsid --w-ld-chr $LDSCORES_DIR/UKBB.EUR.rsid --pop-prev 0.10 --samp-prev 0.344 --out RA_h2_European

# Perform LD score regression for Asian population
ldsc.py --h2 RA_Asian.sumstats.gz --ref-ld-chr $LDSCORES_DIR/UKBB.EAS.rsid --w-ld-chr $LDSCORES_DIR/UKBB.EAS.rsid --pop-prev 0.10 --samp-prev 0.344 --out RA_h2_Asian

# Print out the results for comparison
echo "Heritability estimates:"
echo "European population: $(grep "Heritability" RA_h2_European.log)"
echo "Asian population: $(grep "Heritability" RA_h2_Asian.log)"
