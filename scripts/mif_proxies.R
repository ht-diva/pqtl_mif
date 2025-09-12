

path_base <- "/scratch/dariush.ghasemi/projects/SNPS_IN_LD/output_all/snps_ld_in_meta/"
path_proxies_mif1 <- glue(path_base, "22_24266831_C_G_seq.5356.2_chr22_24241101_24400360_snps_ld_in_pwas_mlog10p.txt")
path_proxies_mif2 <- glue(path_base, "22_24266867_A_G_seq.8221.19_chr22_24234172_24401503_snps_ld_in_pwas_mlog10p.txt")
path_proxies_liftover <- "/scratch/dariush.ghasemi/projects/pqtl_mif/inputs/mif_proxies_liftover.txt"

proxies_headers <- c("SNP", "MLOG10P")
proxies_mif1 <- fread(path_proxies_mif1, col.names = proxies_headers)
proxies_mif2 <- fread(path_proxies_mif2, col.names = proxies_headers)
proxies_liftover <- fread(path_proxies_liftover)

proxies_liftover <- proxies_liftover %>%
  mutate(
    pos_hg37 = gsub("-", ":", hg19),
    pos_hg38 = gsub("chr22-", "chr22:", hg38)
    ) %>%
  select(pos_hg37, pos_hg38)


rbind(
  proxies_mif1 %>% select(SNP),
  proxies_mif2 %>% select(SNP)
  ) %>%
  unique()
