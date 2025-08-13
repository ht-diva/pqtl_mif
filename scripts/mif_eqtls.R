
library(gtexr)

# the largest locus in MIF region 
# build 37 -> 22:23942068-26128449, 
# build 38 -> 22:23599881-25732482

# Can be a GTEx specific ID (e.g. "Whole_Blood") 
# to see valid values or an Ontology ID
get_tissue_site_detail()

# query region in GTEx
res_query <- get_significant_single_tissue_eqtls_by_location(
  tissueSiteDetailId = "Liver",
  start = 23599881,
  end = 25732482,
  chromosome = "chr22",
  datasetId = "gtex_v10",
  .return_raw = FALSE
)

# index and conditional variants at 5 MIF loci
mif_snps <- c(
  "chr22_23924644_C_G_b38",
  "chr22_23924680_A_G_b38",
  "chr22_23992755_A_C_b38",
  "chr22_23998746_A_G_b38",
  "chr22_24225573_A_G_b38",
  "chr22_23995877_G_T_b38",
  "chr22_23970015_A_G_b38",
  "chr22_24235780_C_T_b38",
  "chr22_24247481_C_G_b38",
  "chr22_24499755_A_G_b38"
)

gene_shapes <- rep(c(15,16,17, 18), 8)

res_query %>%
  dplyr::filter(
    #variantId %in% mif_snps
    ) %>%
  dplyr::select(
    chromosome, pos, snpId, tissueSiteDetailId, geneSymbol, nes, pValue
  ) %>% count(geneSymbol) %>% print(n=Inf)
  
# PheWeb plot
res_query %>%
  dplyr::mutate(
    MLOG10P = -log10(pValue),
    # rename gene symbols
    Gene = str_replace(geneSymbol, "ENSG00000224205", "GSTT3P"),
    Gene = str_replace(Gene, "ENSG00000225282", "KLHL5P1"),
    Gene = str_replace(Gene, "ENSG00000228039", "GSTT2B-DDTL"),
    Gene = str_replace(Gene, "ENSG00000231466", "PHF10P2"),
    Gene = str_replace(Gene, "ENSG00000272578", "GUSBP11"),
    Gene = str_replace(Gene, "ENSG00000284128", "BCRP3"),
    Gene = str_replace(Gene, "ENSG00000284233", "GGT5"),
    #Gene = str_replace(Gene, "", ""),
    ) %>% #count(Gene) %>% print(n=Inf)
  #dplyr::filter(pValue < 5e-8) %>%
  ggplot(aes(x = pos, y = MLOG10P, shape = Gene, color = Gene)) +
  geom_point(alpha = 0.9) +
  scale_shape_manual(values= gene_shapes) +
  scale_x_continuous(
    breaks = c(seq(23500000, 25600000, 100000)), 
    labels = c(seq(23.5, 25.6, .1))
    )+
  scale_y_continuous(breaks = seq(0, 80, 10), limits = c(0, 80))+
  labs(
    x = "Genomic position of variant (Kbp)",
    #fill = "Expressed Gene", 
    )+
  #guides(shape = guide_legend(title = "Expressed Gene")) +
  theme_classic() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.85, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    #legend.margin = margin(6, 6, 6, 6),
    #axis.title.x = element_blank()
    axis.ticks.length = unit(2, "mm")
  )


ggsave(
  filename = "13-Aug-25_query_of_mif_in_gtex_v10.png", 
  height = 7.5, width = 13, dpi = 300
)


#===================---==========#
# ----  query SNPs in GTEx  ----
#======================--========#

path_query <- "/scratch/dariush.ghasemi/projects/pqtl_mif/scripts/12-Aug-25_summary.tsv"

# read results for 10 interrogated variants
mif_query <- fread(path_query)

mif_query %>%
  dplyr::mutate(Gene = str_replace(Gene, "ENSG00000224205", "GSTT3P")) %>%
  count(Variant_id)



