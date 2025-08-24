#!/bin/bash

#SBATCH --job-name plink
#SBATCH --output %j_ld_matrix.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1
#SBATCH --mem 8G
#SBATCH --time 1-00:00:00

# imputed genotype of INTERVAL in bed/bim/fam
genotype="/exchange/healthds/pQTL/INTERVAL/Genetic_QC_files/bed/qc_recoded_harmonised/impute_recoded_selected_sample_filter_hq_var_new_id_alleles_"

# load plink
source /exchange/healthds/singularity_functions


# MIF signal
#seqid="seq.8221.19"
#locus="22_24234172_24401503"

# GGT1 signal
seqid="seq.21548.20"
locus="22_23942068_26128449"


# extract params
chr=$(echo $locus | cut -d_ -f1)
beg=$(echo $locus | cut -d_ -f2)
end=$(echo $locus | cut -d_ -f3)

# create output directory
mkdir -p results/bed
mkdir -p results/ld
mkdir -p results/ldmat

# subset the genotype
plink --bfile ${genotype}${chr} --chr ${chr} --from-bp ${beg} --to-bp ${end} --make-bed --out results/bed/${locus} --threads 8  --memory 7500

# compute LD in long format
#plink --bfile results/${locus} --r2 --ld-window 99999 --ld-window-kb 500 --ld-window-r2 0 --out results/ld/${locus} --threads 8  --memory 16000

# extract SNP positions from .bim file
awk '{print $2, $4}' results/bed/${locus}.bim > results/ldmat/${locus}.snps

# compute LD --> '' flag cannot be used with the --r2 matrix
plink --bfile results/bed/${locus} --r2 square --out results/ldmat/${locus} --threads 8  --memory 7500
