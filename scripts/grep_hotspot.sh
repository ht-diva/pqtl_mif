#!/bin/bash

#SBATCH --job-name plink
#SBATCH --output %j_ld_locuszoomr.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1
#SBATCH --mem 8G
#SBATCH --time 1-00:00:00

# imputed genotype of INTERVAL in bed/bim/fam
genotype="/exchange/healthds/pQTL/INTERVAL/Genetic_QC_files/bed/qc_recoded_harmonised/impute_recoded_selected_sample_filter_hq_var_new_id_alleles_"

# load plink
source /exchange/healthds/singularity_functions


# Hotspot CHR2
seqid="seq.21548.20"
locus="2_3600000_3800000"
snp="2:3691548:A:G"

# extract params
chr=$(echo $locus | cut -d_ -f1)
beg=$(echo $locus | cut -d_ -f2)
end=$(echo $locus | cut -d_ -f3)

# create output directory
mkdir -p results/bed
mkdir -p results/ld
mkdir -p results/ldmat


# subset the genotype
# plink --bfile ${genotype}${chr} \
#       --chr ${chr} --from-bp ${beg} --to-bp ${end} \
#       --make-bed \
#       --out results/bed/${locus} \
#       --threads 8 \
#       --memory 7500

#results/${locus} \
# compute LD in long format
plink --bfile ${genotype}${chr} \
      --keep-allele-order \
      --ld-snp $snp \
      --chr $chr  --from-bp $beg  --to-bp $end \
      --r2 \
      --ld-window 99999 \
      --ld-window-kb 500 \
      --ld-window-r2 0 \
      --out results/ld/${locus} \
      --threads 8 \
      --memory 7500
