

path_map <- "/exchange/healthds/pQTL/Reference_datasets_for_QC_proteomics/Cis_trans_mapping/somascan_tss_ncbi_grch37_ensembl_version_20241216.txt"
# inputs (locus_breaker or LB results)
path_freez <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/Locus_breaker_cojo_frozen_version_1812024/"
path_lb_cistrans <- "Archived/mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann.csv"
path_cojo <- "16-Dec-24_collected_independent_snps.csv"


# read data files
lb_cistrans <- fread(paste0(path_freez, path_lb_cistrans))
cojo <- read.csv(paste0(path_freez, path_cojo))
map <- data.table::fread(path_map)



mif_pqtls <- c('seq.19230.12', 'seq.19237.17', 'seq.21548.20', 'seq.4878.3', 'seq.5356.2', 'seq.8221.19')

# Find TSS and gene name for pQTLs in the broader MIF region
mif_genes <- map %>% 
  dplyr::filter(target %in% mif_pqtls) %>%
  dplyr::select(
    target, symbol, chr_TSS = chromosome, TSS, #Target_Name, Target_Full_Name
    )

# loci or pQTLs fell into broader MIF region
mif_loci <- fread(out_lb_epitop_cojo ) %>% #lb_cistrans %>%
  dplyr::filter(
    chr == 22, POS >= 23942068, POS <= 26128449
  ) %>%
  dplyr::select(
    chr:SNPID, MLOG10P, phenotype_id, cis_or_trans
  )

mif_cojo <- cojo %>% 
  dplyr::filter(Chr == 22, study_id %in% mif_pqtls) %>%
  group_by(study_id, locus) %>%
  summarize(
    cojo_snps = paste0(unique(SNP), collapse = "; ")
  ) %>%
  ungroup()

mif_loci %>%
  left_join(
    mif_genes,
    join_by(phenotype_id == target)
  ) %>%
  dplyr::mutate(
    distance = abs(POS - TSS) # gap between pQTL and TSS
  ) %>%
  left_join(
    mif_cojo,
    join_by(phenotype_id == study_id)
  ) %>%
  dplyr::select(
    cis_trans = cis_or_trans,
    seqid = phenotype_id,
    symbol,
    chr:end,
    SNPID,
    MLOG10P,
    chr_TSS,
    TSS,
    distance,
    cojo_snps
  )
