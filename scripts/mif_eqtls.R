library(tidyverse)
library(gtexr)

# inputs
path_query <- "/scratch/dariush.ghasemi/projects/pqtl_mif/results/12-Aug-25_summary.tsv"

# outputs
tbl_region_query <- "/scratch/dariush.ghasemi/projects/pqtl_mif/results/22-Aug-25_region_eqtls_interrogation_gtex_v10.csv"
scatter_eqtls <- "13-Aug-25_query_of_mif_in_gtex_v10.png"
scatter_pdf <- "23-Aug-25_query_of_mif_in_gtex_v10.pdf"
heatmap_eqtls <- "05-Sep-25_query_of_mif_proxies_in_gtex_v8.png"

#--------------------------------# 
# the largest locus in MIF region 
# build 37 -> 22:23942068-26128449
# build 38 -> 22:23599881-25732482

# Gene region: GRCh38 <-> GRCh37
#  MIF: 23894383..23895223 <-> 24236570..24237410
#  DDT: 23971370..23980504 <-> 24313559..24322695
#--------------------------------# 

# Can be a GTEx specific ID (e.g. "Whole_Blood") 
# to see valid values or an Ontology ID
get_tissue_site_detail()

# query region in GTEx
res_query_v8 <- get_significant_single_tissue_eqtls_by_location(
  tissueSiteDetailId = "Liver",
  start = 23599881,
  end = 25732482,
  chromosome = "chr22",
  datasetId = "gtex_v8",
  .return_raw = FALSE
)

annot_eqtl <- data.frame(
         snpId = c("rs5760119", "rs5760120", "rs6004011", "rs4822466"),
         locus = c("MIF index", "MIF index", "GSTT1 conditional", "DDT conditional")
         )

#--------------------------------#
# Visualize query results for mif index/cojo snps
res_query_v8 %>%
  dplyr::select(
    chromosome, pos, snpId, tissueSiteDetailId, geneSymbol, nes, pValue
  ) %>%
  left_join(annot_eqtl, join_by(snpId)) %>%
  #count(geneSymbol) %>% print(n=Inf)
  mutate(
    chr_pos = str_c(chromosome, ":", pos),
    my_snp  = paste0(snpId, ", ", locus),
    my_snp  = str_remove(my_snp, ", NA"),
    snp_annot = fct_reorder(my_snp, pos)
    ) %>%
  dplyr::filter(
    chr_pos %in% proxies_liftover$pos_hg38, # of 39 proxies, v8 had 34 eQTLs; v10 had 35.
    pValue < 5e-8
    ) %>%
  ggplot(aes(geneSymbol, snp_annot, fill = nes)) +
  geom_tile(color = "grey30") +
  scale_y_discrete(position = "right")+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000")+ #coord_fixed()
  guides(fill = guide_colourbar(position = "bottom", title = "Normalized\neffect size"))+
  theme_bw()+
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12, face = 3)
  )

# save heatmap
ggsave(filename = heatmap_eqtls, height = 7.5, width = 13, dpi = 150)


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



