#!/bin/bash

#SBATCH --job-name grep_snps
#SBATCH --output %j_extract_cs.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1
#SBATCH --mem 2G
#SBATCH --time 1-00:00:00

# imputed genotype of INTERVAL in bed/bim/fam
genotype="/exchange/healthds/pQTL/INTERVAL/Genetic_QC_files/bed/qc_recoded_harmonised/impute_recoded_selected_sample_filter_hq_var_new_id_alleles_"

# load plink
source /exchange/healthds/singularity_functions

chrom=22

# create output directory
#mkdir -p results/dosage


# compute LD --> '' flag cannot be used with the --r2 matrix
plink --bfile $genotype${chrom} \
      --keep-allele-order \
      --extract ../inputs/mif_eqtls_credible.snps \
      --recode A  \
      --out ../results/dosage/mif_eqtls_cs \
      --threads 8  \
      --memory 7500

