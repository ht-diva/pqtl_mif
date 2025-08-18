
# inputs
path_query <- "/scratch/dariush.ghasemi/projects/pqtl_mif/results/gwas_catalog_22_23824644_to_24024644.tsv"


query_res <- fread(path_query, sep = "\t")

query_res %>%
  dplyr::select(
    CHR_ID,
    CHR_POS,
    SNPS,
    reported_gene = "REPORTED GENE(S)",
    MAPPED_GENE,
    disease_trait = "DISEASE/TRAIT",
    MAPPED_TRAIT,
    PVALUE_MLOG,
    effect_size = "OR or BETA"
  ) %>% 
  dplyr::filter(PVALUE_MLOG > 7.3) %>% 
  dplyr::mutate(
    # gene = str_replace_all(MAPPED_GENE, "MIF.*$", "MIF"),
    # gene = str_replace_all(gene, "GSTT3P.*$", "GSTT3P"),
    # gene = str_replace_all(gene, "GSTT2.*$", "GSTT2"),
    # gene = str_replace_all(gene, "KLHL5P1.*$", "GSTT2"),
    # gene = str_replace_all(gene, ".*DDTL.*$", "DDT"),
    # gene = str_replace_all(gene, ".*SLC2A11", "SLC2A11"),
    # gene = str_replace_all(gene, ".* DERL3|DERL3.*$", "SMARCB1"),
    # gene = str_replace_all(gene, ".* CABIN1|CABIN1.*$", "GSTT4"),
    trait = str_remove_all(MAPPED_TRAIT, "measurement| serum | serum|serum | in blood| blood serum | protein level ratio"),
    trait = str_replace_all(trait, "coiled-coil-helix-coiled-coil-helix domain-containing protein 10.*$", "CHCHD10"),
    trait = str_replace_all(trait, "scavenger receptor cysteine-rich domain-containing group B protein", "SRCR-SF"),
    trait = str_replace_all(trait, "IgG .*$", "IgG related"),
    trait = str_replace_all(trait, "high density lipoprotein .*$", "HDL-C"),
    trait = str_replace_all(trait, "gamma-glutamyl transferase .*$", "GGT"),
    trait = str_replace_all(trait, "alkaline phosphatase .*$", "ALP"),
    trait = str_replace_all(trait, "aspartate aminotransferase .*$", "AST"),
    trait = str_replace_all(trait, "alanine aminotransferase .*$", "ALT"),
    trait = ifelse(disease_trait == "Glutathione S-transferase theta-1 levels", "Glutathione S-transferase theta-1 levels", trait),
    trait = fct_reorder(trait, -PVALUE_MLOG), 
    ) %>%
  ggplot(aes(x = trait, y = PVALUE_MLOG, color = MAPPED_GENE)) + #, shape = gene
  geom_point(size = 3) +
  #scale_shape_manual(values = c(19, 18, 19, 17, 19, 17, 17, 19))+
  #scale_color_brewer(palette = "Paired")+
  #scale_y_discrete(position = "right") +
  #paletteer::scale_colour_paletteer_d("vapoRwave::hotlineBling")+
  ggtitle("Query of MIF region in GWAS Catalog") +
  theme_classic()+
  theme(
    legend.position = c(.9, .6),
    axis.text.x = element_text(size=8, face = 1, hjust = 0, vjust = 0.5, angle = -90),
    axis.title.x = element_blank(),
    plot.margin = margin(r = 2, b = 2, t = 2, l = 2, unit = "mm")
  )

ggsave(
  filename = "04-Aug-25_query_of_mif_in_gc.png", 
  height = 5.5, width = 8, dpi = 300
  )

