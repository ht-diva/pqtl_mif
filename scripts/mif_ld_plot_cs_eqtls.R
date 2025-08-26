
library(tidyverse)
library(data.table)
library(glue)
library(corrplot)
library(pheatmap)

#--------------------------------# 
# inputs
path_base <- "/scratch/dariush.ghasemi/projects/pqtl_mif/"
path_dosage_mif_index <- glue(path_base, "results/dosage/mif_cojo.raw")
path_dosage_eqtls_cs  <- glue(path_base, "results/dosage/mif_eqtls_cs.raw")

# outputs
heatmap_sc_eqtls <- "25-Aug-25_LD_index_cojo_vs_eqtls_cs.pdf"


#--------------------------------# 
# Gene used for annotation of heatmap
gstt1_snps <- c("22:24334948:A:C", "22:24338071:G:T", "22:24340940:A:G", "22:24398868:A:C")
ggt5_snps  <- c("22:24621541:A:G", "22:24631748:C:T", "22:24643449:C:G", "22:24895723:A:G")
ddt_snps   <- c("22:24312204:A:G", "22:24340940:A:G")
f10_snps   <- "22:24997846:C:T"

#--------------------------------# 
# Read genotype file
geno_index <- fread(path_dosage_mif_index)
geno_eqtls <- fread(path_dosage_eqtls_cs)


#=============================#
# -----    MIF index     -----
#=============================#

# calculate LD between index variants of the 6 loci in broader MIF region 

# remove effective allele from column names
geno <- geno_index %>% 
  inner_join(geno_eqtls) %>% # join dosage of GTEx credible set eQTLs
  rename_with(~ str_remove(., "_[ATCG]+$"), matches("_[ATCG]+$"))

# store variants list for the match
snps <- names(geno)[7:ncol(geno)]

# convert geneotype to matrix
X <- as.matrix(geno[,7:ncol(geno)]) * 1.0

# compute LD (r2)
R <- cor(X, use = "complete.obs")  # ---> Find a solution for NAs in genotype and cor matrix

# corrplot.mixed(R)
# 
# corrplot(
#   R, 
#   method = 'number', 
#   diag = FALSE,
#   #type = 'upper',
#   col = COL2('PRGn') #COL2('BrBG')
# ) # colorful number


#=============================#
# ----- heatmap annotation ----
#=============================#

# merge cs variants positions (b37&b38) to their corresponding eQTL gene 
annot_heatmap <- cs_snps_b37 %>% 
  #dplyr::filter(str_detect(chr_pos_b37, "22:24266831|22:24266867")) %>%
  left_join(res_cs_snps, join_by(variantId)) %>% # join cs variants with eGenes
  dplyr::mutate(
    chr_pos = str_remove(chr_pos_b37, ":[ATCG]+:[ATCG]+$"),
    gene = str_replace_all(gencodeId, "ENSG00000099984.11", "GSTT2"),
    gene = str_replace_all(gene, "ENSG00000128262.8",  "POM121L9P"),
    gene = str_replace_all(gene, "ENSG00000189269.12", "DRICH1"),
    gene = str_replace_all(gene, "ENSG00000218537.1",  "MIF-AS1")
  ) %>% 
  dplyr::select(chr_pos, gene) %>% # take variants positions in b37 and eGenes
  right_join(
    # transform correlation matrix to df and extract chromosomal positions for join
    data.frame(chr_pos_geno = colnames(R)) %>%
      dplyr::mutate(chr_pos = str_remove(chr_pos_geno, ":[ATCG]+:[ATCG]+$")),
    join_by(chr_pos)
  ) %>%
  dplyr::mutate(  # define eGene or annotated gene and type of variant (index/cojo/eQTL cs)
    gene = ifelse(chr_pos_geno %in% gstt1_snps, "GSTT1", gene),
    gene = ifelse(chr_pos_geno %in% ggt5_snps, "GGT5", gene),
    gene = ifelse(chr_pos_geno %in% ddt_snps, "DDT", gene),
    gene = ifelse(chr_pos_geno %in% f10_snps, "F10", gene),
    variant = ifelse(chr_pos_geno %in% c(gstt1_snps, ggt5_snps, ddt_snps, f10_snps), "index/cojo", "eQTL cs"),
    variant = ifelse(str_detect(chr_pos, "22:24266831|22:24266867"), "index + eQTL cs", variant)
  ) %>%
  column_to_rownames("chr_pos_geno") %>% # to match with row-names in correlation matrix
  dplyr::select(variant, gene)


#=============================#
# -----   heatmap  plot   ----
#=============================#

# draw heatmap of LD
pdf(heatmap_cs_eqtls, width = 18, height = 15) #, units = "cm", res = 300

pheatmap(R, annotation_row = annot_heatmap)

dev.off()

