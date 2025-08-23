library(tidyverse)
library(gtexr)

# inputs
path_query <- "/scratch/dariush.ghasemi/projects/pqtl_mif/results/12-Aug-25_summary.tsv"

# outputs
tbl_region_query <- "/scratch/dariush.ghasemi/projects/pqtl_mif/results/22-Aug-25_region_eqtls_interrogation_gtex_v10.csv"
scatter_eqtls <- "13-Aug-25_query_of_mif_in_gtex_v10.png"
scatter_pdf <- "23-Aug-25_query_of_mif_in_gtex_v10.pdf"

#--------------------------------# 
# the largest locus in MIF region 
# build 37 -> 22:23942068-26128449
# build 38 -> 22:23599881-25732482

# Gene region: GRCh38 <-> GRCh37
#  MIF: 23894383..23895223 <-> 24236570..24237410
#  DDT: 23971370..23980504 <-> 24313559..24322695
#--------------------------------# 

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
#--------------------------------# 

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


res_query %>%
  dplyr::filter(
    #variantId %in% mif_snps
    ) %>%
  dplyr::select(
    chromosome, pos, snpId, tissueSiteDetailId, geneSymbol, nes, pValue
  ) %>% count(geneSymbol) %>% print(n=Inf)

#--------------------------------#
# QC query results
res_qced <- res_query %>%
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
    )# %>% count(Gene) %>% print(n=Inf)


# save results for Adam
write.csv(res_qced, tbl_region_query, row.names = F, quote = F)

# find cis-eQTLs around DDT-MIF genes
res_qced %>%
  dplyr::filter(
    pValue < 5e-8,
    pos > 23894383 & pos < 23980504
    ) %>%
  count(variantId)



#--------------------------------# 
# PheWeb plot

gene_shapes <- rep(c(15,16,17, 18), 8)

plt_eqtls <- res_qced %>%
  dplyr::filter(
    pValue < 5e-8,
    pos > 23894383 & pos < 23980504
    ) %>%
  ggplot(aes(x = pos, y = MLOG10P, shape = Gene, color = Gene)) +
  geom_point(alpha = 0.9) +
  scale_shape_manual(values= gene_shapes) +
  scale_x_continuous(
    breaks = c(seq(23890000, 23980000, 10000)),
    # breaks = c(seq(23500000, 25600000, 100000)),
    labels = c(seq(23.89, 23.98, .01)),
    # labels = c(seq(23.5, 25.6, .1))
    limits = c(23890000, 23980000)
    )+
  scale_y_continuous(breaks = seq(0, 80, 10), limits = c(0, 80))+
  labs(
    x = "Genomic position of variant (Mb)",
    #fill = "Expressed Gene", 
    )+
  #guides(shape = guide_legend(title = "Expressed Gene")) +
  theme_classic() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.22, .98),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.background = element_blank(),
    #legend.margin = margin(6, 6, 6, 6),
    #axis.title.x = element_blank()
    axis.ticks.length = unit(2, "mm")
  )


ggsave(filename = scatter_eqtls, height = 7.5, width = 13, dpi = 300)


#--------------------------------# 
#BiocManager::install("EnsDb.Hsapiens.v86")
library(locuszoomr)
library(EnsDb.Hsapiens.v86) # GRCh38
library(EnsDb.Hsapiens.v75) # GRCh37

data("SLE_gwas_sub") # built-in mini dataset


loc <- locus(
  data = SLE_gwas_sub,
  gene = 'AP000350.6', 
  fix_window = 1.1e5,
  ens_db = "EnsDb.Hsapiens.v86"
  # ens_db = "EnsDb.Hsapiens.v75"
  )

# join scatter plot + gene annotation
pdf(scatter_pdf, width = 8, height = 5.5)

# set up layered plot; store old par() settings
oldpar <- set_layers(1)

plt_eqtls
genetracks(loc, highlight = c("MIF", "DDT"))

par(oldpar)  # revert par() settings
dev.off()


#===================---==========#
# ----  query SNPs in GTEx  ----
#======================--========#


# read results for 10 interrogated variants
mif_query <- fread(path_query)

mif_query %>%
  dplyr::mutate(Gene = str_replace(Gene, "ENSG00000224205", "GSTT3P")) %>%
  count(Variant_id)



