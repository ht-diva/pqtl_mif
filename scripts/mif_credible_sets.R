
library(gtexr)

# define involved genes in MIF region
mif_genes <- c(
  "ENSG00000189269.12", # DRICH1
  "ENSG00000099958.15", # DERL3
  "ENSG00000099953.10", # MMP11
  "ENSG00000240972.2",  # MIF
  "ENSG00000218537.1",  # MIF-AS1
  "ENSG00000133460.20", # SLC2A11
  "ENSG00000099977.15", # DDT
  "ENSG00000099974.8",  # DDTL
  "ENSG00000099994.11", # SUSD2
  "ENSG00000099998.19", # GGT5
  "ENSG00000099984.11", # GSTT2
  "ENSG00000133433.11", # GSTT2B
  "ENSG00000276950.6",  # GSTT4
  "ENSG00000099956.20", # SMARCB1
  "ENSG00000099991.18", # CABIN1
  "ENSG00000128262.8",  # POM121L9P
  "ENSG00000225098.1",  # BCRP1
  "ENSG00000128218.8",  # VPREB3
  "ENSG00000100014.20", # SPECC1L
  "ENSG00000178803.12", # ADORA2A-AS1
  "ENSG00000100031.19"  # GGT1
)

# search by gene
res_indep_snps <- get_independent_eqtl(
  gencodeIds = mif_genes, 
  tissueSiteDetailIds = "Liver",
  datasetId = "gtex_v8",
  .return_raw = FALSE,
  itemsPerPage = 100000
  )

res_indep_snps %>% #dim()
  dplyr::filter() %>%
  pull(snpId)



#================================#
#---- credible sets in liver ---- 
#================================#

# optionally filter for a single variant and/or one or more tissues
res_cs_snps <- get_fine_mapping(
  gencodeIds = mif_genes,
  #variantId = "chr1_153228363_A_G_b38",
  tissueSiteDetailIds = "Liver"
)


res_cs_snps %>% dplyr::filter(variantId %in% mif_snps)
